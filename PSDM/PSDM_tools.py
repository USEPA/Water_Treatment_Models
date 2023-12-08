# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 10:07:16 2020

PSDM Tools 
----------

This file contains tools and associated functions geared at aiding in the 
analysis of GAC adsorption data. 

Functions Contained:
    isotherm_fit()
    predict_full_scale()
    specific_throughput()


@author: Jonathan Burkhardt

EPA Disclaimer
==============
The United States Environmental Protection Agency (EPA) GitHub project code is 
provided on an "as is" basis and the user assumes responsibility for its use. 
EPA has relinquished control of the information and no longer has 
responsibility to protect the integrity , confidentiality, or availability of 
the information. Any reference to specific commercial products, processes, or 
services by service mark, trademark, manufacturer, or otherwise, does not 
constitute or imply their endorsement, recomendation or favoring by EPA. The 
EPA seal and logo shall not be used in any manner to imply endorsement of any 
commercial product or activity by EPA or the United States Government.

By submitting a pull request, you make an agreement with EPA that you will not 
submit a claim of compensation for services rendered to EPA or any other 
federal agency. Further, you agree not to charge the time you spend developing 
software code related to this project to any federal grant or cooperative 
agreement.
"""
import pandas as pd
import warnings
warnings.simplefilter("ignore")
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import multiprocessing as mp

import PSDM
from PSDM_functions import find_minimum_df
import time as ti

def specific_throughput(column_specs, filter_pfas, k_data, c0, ct,
                        compound_data,
                        ebct_range=np.arange(5,46,5),
                        chem_type=['halogenated alkenes'],
                        wat_type=['Organic Free'],
                        nr=14, nz=18):
    '''
    

    Parameters
    ----------
    column_specs : pandas dataframe
        Dataframe with column specifications.
    filter_pfas : list of strings
        List of compounds to model.
    k_data : pandas dataframe
        Dataframe that contains K & 1/n values to be used in the simulation.
        Must provide a k_data structure that has all the compounds listed in 
        filter_pfas. Can contain more species, but must contain all requested 
        by filter_pfas.
    c0 : float
        Initial concentration of chemicals in influent. Units: ng/L
        Assumes all chemicals have the same initial concentration.
        Does not allow for variable influent concentrations.
    ct : float
        Threshold concentration for when a bed is removed. Units: ng/L
        Defines carbon usage rate. Must be less than c0.
    compound_data : pandas dataframe
        Dataframe that contains physical parameters associated with chemical
        species included in simulation.
    ebct_range : list/iterable, optional
        Values of empty bed contact time to consider. 
        The default is np.arange(5,46,5).
    chem_type : list/string, optional
        Type of chemical species to model. The default is ['halogenated alkenes'].
        Related to fouling parameters
    wat_type : list/string, optional
        Type of water to model. The default is ['Organic Free'].
        Related to fouling parameters
    nr : int, optional
        Number of radial collocation points. The default is 14.
    nz : int, optional
        Number of axial collocation points. The default is 18.

    Returns
    -------
    compound_store : TYPE
        DESCRIPTION.

    '''
    
    orig_ebct = column_specs['L'] * np.pi * (column_specs['diam']**2)/\
                (4. * column_specs['flrt'])
    orig_flrt = column_specs['flrt'] * 1
    orig_L = column_specs['L'] * 1.
    orig_diam = column_specs['diam'] * 1.
    orig_wt = column_specs['wt'] * 1.
    
    types = [column_specs['carbon'], column_specs['influentID']]
    multi_idx = pd.MultiIndex.from_tuples([(typ, comp)
                                           for typ in types 
                                           for comp in filter_pfas],
                                          names=['type','compound'])
    idx = [0, column_specs['duration']]
    raw_data = pd.DataFrame(c0, columns=multi_idx, index=idx)

    #Initiate storage dictionary (returned object)
    compound_store = {}
   
    for comp in filter_pfas:
        print(comp)
        ebct_store = []
        for ebct in ebct_range:
            ratio = orig_ebct / ebct
            #rescale flow rate of system to desired EBCT value
            if ebct <= orig_ebct:
                column_specs['flrt'] = ratio * orig_flrt
                column_specs['diam'] = orig_diam * 1.
                column_specs['L'] = orig_L * 1.
                column_specs['wt'] = orig_wt
            else:
                #Resize bed, not adjust flowrate, should allows shorter durations
                column_specs['flrt'] = orig_flrt * 1.
                r_ratio = ratio**(1/3.)
                column_specs['diam'] = orig_diam * r_ratio
                column_specs['L'] = orig_L * r_ratio
                column_specs['wt'] = orig_wt * ratio # not modified ratio
                
            
            #need to rework this to support this step...
            column = PSDM.PSDM(column_specs, compound_data, raw_data,\
                            nz=nz, nr=nr, chem_type=chem_type,\
                            water_type=wat_type, k_data=k_data,\
                            xn_range=[k_data[comp]['1/n']],
                            test_range=[k_data[comp]['K']],
                            optimize=False)
            
            print(column.ebct)
            _, _, _, _, results = column.run_psdm_kfit(comp)
            
            treat_days = results[results < ct].dropna().index[-1]
            spec_throughput = (column.flrt/1e6 * PSDM.min_per_day * \
                              treat_days) / (column.wt/1e3)
            
            ebct_store.append(spec_throughput)
        
        compound_store[comp] = ebct_store
    
    return compound_store

def predict_full_scale(PSDM_obj, filter_pfas, target_conc, \
                       total_beds, beds_per_cycle, plot=True):
    '''
    

    Parameters
    ----------
    PSDM_obj : PSDM class object
        Column information created from PSDM.PSDM()
        Must have 'k_data=' supplied on object creation. Should only use 
        user-supplied k_values (or initial estimates will be used)
    filter_pfas : list of strings
        Example: ['compound1', 'compound2',...]
        List of compounds to model and use to esablish the target_conc.
        if only a single compound is needed : ['compound'] must be supplied
    target_conc : float
        The target concentration for cummulative modeled effluent from 
        filter_pfas. Units are in ng/L (ppt).
    total_beds : INT
        Number of beds in rotation.
    beds_per_cycle : INT
        Number of beds rotated in/out for each cycle.
    plot : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    best_val : float
        Number of days per bed rotation interval. This highlights how many 
        days between bed replacements. 
        
        Example: best_val = 100 (days), 
                 for 8 total beds and 2 beds per cycle
                 means 2 beds are cycled in every 100 days, for total life
                 of any 1 bed of 400 days (8/2 * 100 days)
        
    best_cycle : interpolating function
        Blended effluent concentration for best_val case.
        y = concentration (ng/L)
        x = time (days)
        
    Plots and Excel files are also generated

    '''
    init_optflag = PSDM_obj.optimize_flag
    init_testrange = PSDM_obj.test_range
    init_xnrange = PSDM_obj.xn_range
    init_compounds = PSDM_obj.compounds
    init_xdata = PSDM_obj.xdata
    init_xn = PSDM_obj.xn

    
    PSDM_obj.optimize_flag = False
    PSDM_obj.compounds = filter_pfas
    
    idx = PSDM_obj.data_df.index.values
    if np.max(idx) < PSDM_obj.max_days:
        idx[-1] = PSDM_obj.max_days
        PSDM_obj.data_df.set_index(idx)
        #this assumes that the concentrations in the dataframe are just 
        #averages, so no time variability is impacted
    
    PSDM_obj.xdata = PSDM_obj.data_df.index.values
    
    time_idx = np.arange(PSDM_obj.max_days+1)
    data_df = pd.DataFrame(columns=filter_pfas, index=time_idx)
    
    for comp in filter_pfas:
        PSDM_obj.test_range = np.array([PSDM_obj.k_data[comp]['K']])
        PSDM_obj.xn_range = np.array([PSDM_obj.k_data[comp]['1/n']])
        PSDM_obj.xn = PSDM_obj.k_data[comp]['1/n']
        comp, k_v, xn_v, ssqs, md = PSDM_obj.run_psdm_kfit(comp)

        md[md<0.] = 0.
        md[md>data_df[comp][0]] = data_df[comp][0]
        if plot:
            plt.plot(md.index.values, md.values, label=comp)
        out_f = interp1d(md.index.values,\
                         md.transpose().values,\
                         fill_value='extrapolate')
        data_df[comp] = out_f(time_idx)[0]
        
    if plot:
        plt.legend(loc='center right')
        plt.ylabel('Concentration (ng/L)')
        plt.xlabel('Time (days)')
        plt.xlim((0,1095)) #limits to 3 years, rather than 3000 days
        #may change to 1000 days
        plt.savefig('full_scale_'+PSDM_obj.carbon+'.png',dpi=300)
        plt.close()
    
    data_df[data_df<0]=0. #resets negatives to zero
    
    writer = pd.ExcelWriter(PSDM_obj.project_name+'_'+PSDM_obj.carbon+'.xlsx')
    data_df.to_excel(writer, 'model_fit')
    writer.save()

    small_rotation = total_beds%beds_per_cycle
    
    if small_rotation == 0:
        #all cycles the same
        weights = int(total_beds/beds_per_cycle) *\
                  [float(beds_per_cycle)/total_beds]
    else:
        weights = [float(small_rotation)/total_beds] + \
                  int(total_beds/beds_per_cycle)*[float(beds_per_cycle)/total_beds]
        # having small rotation be first, is the worst case scenario
        # it means the smallest percentage of beds is new
        # having small rotation last is the best case scenario, not used.
    
    summed = data_df.transpose().sum()
    
    min_cycle = 1 #days
    
    num_cycles = np.ceil(float(total_beds)/beds_per_cycle)
    if num_cycles*beds_per_cycle > total_beds:
        print('number of beds per cycle may result in non-uniform cycling. assuming new beds have fewest numbers')
    
    bed_info = []
    count = 0
    bed_c = 1
    for bed in range(total_beds):
        if count <= beds_per_cycle:
           bed_info.append(bed_c)
        count+=1
        if count == beds_per_cycle:
            count = 0
            bed_c+=1
           
    function = interp1d(summed.index,summed.values,fill_value='extrapolate')
    aa = np.arange(5*PSDM_obj.max_days)
    summed = pd.Series(function(aa),index = aa)
    
    best_cycle = np.zeros(min_cycle)
    best_val = min_cycle*1
    try:
        for i in range(min_cycle,PSDM_obj.max_days,1):
            tmp = np.zeros(i)
            for j in range(max(bed_info)):
                tmp += weights[j]*(summed[(summed.index>=(j*i))&(summed.index<((j+1)*i))].values)
            if tmp.max() <= target_conc:
                best_cycle = tmp*1.
                best_val = i
            else:
                break
    except Exception:
        best_val = PSDM_obj.max_days
        best_cycle = np.zeros(PSDM_obj.max_days)

    #reset object parameters to initial values
    PSDM_obj.optimize_flag = init_optflag 
    PSDM_obj.test_range = init_testrange 
    PSDM_obj.xn_range = init_xnrange 
    PSDM_obj.compounds = init_compounds
    PSDM_obj.xdata = init_xdata
    PSDM_obj.xn = init_xn
    
    return best_val, best_cycle


def isotherm_fit(data, isotherm='freundlich', plot=True, save_plot=False, filename='test'):
    '''
    

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    isotherm : TYPE, optional
        DESCRIPTION. The default is 'freundlich'.
    plot : TYPE, optional
        DESCRIPTION. The default is True.
    save_plot : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    def langmuir(c, K, N):
        '''
        Returns results of Langmuir Isotherm
        
        Parameters
        ----------
        c : array
            array of liquid concentrations.
        K : float
        N : float
            K & N are parameter values
        
        Returns
        -------
        array, of solid phase concentrations

        '''
        return (K*c*N)/(1. + K*c)
    
    def freundlich(c, k, invN):
        '''
        

        Parameters
        ----------
        c : array
            array of liquid concentrations.
        k : float
            Freundlich K parameter.
        invN : float
            1/n parameter 

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        # k, invN = array
        return k * c**invN
    
    def RedlichPeterson(c, A, B, M):
        return (A*c)/(B+c**M)

    
    # Function below, isotherm equations above 

    if 'q' not in data.columns:
        data['q'] = (data['C0']-data['Ce'])/data['mass']

    # print(data) #debug, to remove

    xdata = data['Ce'].values
    ydata = data['q'].values
    
    xdata = xdata[~np.isnan(ydata)]
    ydata = ydata[~np.isnan(ydata)]
    
    xdata_plot = np.linspace(xdata.min(), xdata.max(), 30)
    
    if plot:
        plt.figure(figsize=(8,5))
        ax = plt.gca()
        plt.plot(xdata, ydata, 'x', label='Experimental')
        
    isotherm_available = True
    if isotherm.lower() == 'freundlich':
        print('Freundlich: K * Ce**(1/n)')
        popt, pcov = curve_fit(freundlich, xdata, ydata, 
                               bounds=(0, [np.inf, np.inf]),
                               maxfev=10000) 
        intervals = np.sqrt(np.diag(pcov))
        tmp_data = pd.DataFrame(index=['K','1/n'], columns=['parameter', 'error'])
        tmp_data['parameter'] = popt
        tmp_data['error'] = intervals
        y_plot = freundlich(xdata_plot, popt[0], popt[1])
        y_model = freundlich(xdata, popt[0], popt[1])
        
        plot_title = (popt[0], intervals[0], popt[1], intervals[0])
        title = 'Freundlich\nK: %.3e$\pm$%.3e - 1/n: %.3e$\pm$%.3e' %plot_title
        
    elif isotherm.lower() == 'langmuir':
        print('Langmuir: qm * KL * Ce/(1 + KL*Ce)')
        popt, pcov = curve_fit(langmuir, xdata, ydata, 
                                bounds=(0, [np.inf, np.inf]),
                               maxfev=10000)    
        intervals = np.sqrt(np.diag(pcov))

        tmp_data = pd.DataFrame(index=['KL','qm'], columns=['parameter', 'error'])
        tmp_data['parameter'] = popt
        tmp_data['error'] = intervals

        y_plot = langmuir(xdata_plot, popt[0], popt[1])
        y_model = langmuir(xdata, popt[0], popt[1])
        
        plot_title = (popt[0], intervals[0], popt[1], intervals[0])
        title = 'Langmuir\nK$_{L}:$ %.3e$\pm$%.3e - q$_{m}$: %.3e$\pm$%.3e' %plot_title
    elif isotherm.lower() == 'redlichpeterson':
        popt, pcov = curve_fit(RedlichPeterson, xdata, ydata, 
                               bounds=(0, [np.inf, np.inf, np.inf]),
                               maxfev=10000)    
        intervals = np.sqrt(np.diag(pcov))

        tmp_data = pd.DataFrame(index=['A','B','M'], columns=['parameter', 'error'])
        tmp_data['parameter'] = popt
        tmp_data['error'] = intervals

        y_plot = RedlichPeterson(xdata_plot, popt[0], popt[1], popt[2])
        y_model = RedlichPeterson(xdata, popt[0], popt[1], popt[2])
        
        plot_title = (popt[0], intervals[0], popt[1], intervals[0], popt[2], intervals[2])
        title = 'Redlich Peterson\nA: %.2e$\pm$%.2e - B: %.2e$\pm$%.2e - M: %.2e$\pm$%.2e' %plot_title

    else:
        print('Warning: Isotherm Selection does not match available choices')
        isotherm_available = False
    
    if plot:
        plt.plot(xdata_plot, y_plot , label='Best Fit')

    
    if isotherm_available:
        m = popt.size
        n = xdata.size
        dof = n - m
        
        t = stats.t.ppf(0.975, dof)
        resid = ydata - y_model
        chi2 = np.sum((resid/y_model)**2)
        chi2_red = chi2/dof
        s_err = np.sqrt(np.sum(resid**2)/dof)
        
        ci = t*s_err*np.sqrt(1/n + (xdata_plot-np.mean(xdata))**2/\
                             np.sum((xdata-np.mean(xdata))**2))
        pi = t*s_err*np.sqrt(1+1/n+(xdata_plot-np.mean(xdata))**2/\
                             np.sum((xdata-np.mean(xdata))**2))
        
        if plot:
            ax.fill_between(xdata_plot, y_plot-ci, y_plot+ci, 
                            color = '#b9cfe7', edgecolor = '', 
                            label = '95% Confidence Interval' )
            ax.fill_between(xdata_plot, y_plot-pi, y_plot+pi, 
                            linestyle = '--', color = 'None')
            plt.plot(xdata_plot,y_plot-pi, linestyle = '--', 
                     color = '0.5', label = '95% Prediction Interval')
            plt.plot(xdata_plot,y_plot+pi, linestyle = '--', color = '0.5')
            plt.title(title)
    
    if plot:
        plt.legend()
        plt.xlabel('C$_{e}$: Equilibrium Concentration')
        plt.ylabel('q: Solid Phase Concentration')
        
    if plot and save_plot:
        plt.savefig(filename + '.png') 

    print(tmp_data)
    return tmp_data


# =============================================================================
# NEED TESTING and/or Example file demonstrating use
# =============================================================================
# =============================================================================
# alternate full scale for ds
# 
# =============================================================================
def predict_full_scale_ds(PSDM_obj, filter_pfas, target_conc, \
                       total_beds, beds_per_cycle, plot=True):
    '''
    

    Parameters
    ----------
    filter_pfas : TYPE
        DESCRIPTION.
    target_conc : TYPE
        DESCRIPTION.
    total_beds : TYPE
        DESCRIPTION.
    beds_per_cycle : TYPE
        DESCRIPTION.
    plot : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    best_val : TYPE
        DESCRIPTION.
    best_cycle : TYPE
        DESCRIPTION.

    '''
    
    opt_flg = PSDM_obj.optimize_flag
    PSDM_obj.optimize_flag = False
    
    compounds = PSDM_obj.compounds
    PSDM_obj.compounds = filter_pfas
    
    idx = PSDM_obj.data_df.index.values
    if np.max(idx) < PSDM_obj.max_days:
        idx[-1] = PSDM_obj.max_days
        PSDM_obj.data_df.set_index(idx)
        #this assumes that the concentrations in the dataframe are just 
        #averages, so no time variability is impacted
    
    test_range = PSDM_obj.test_range
    xdata = PSDM_obj.xdata 
    
    PSDM_obj.xdata = PSDM_obj.data_df.index.values
    
    time_idx = np.arange(PSDM_obj.max_days+1)
    data_df = pd.DataFrame(columns=filter_pfas, index=time_idx)
        
    for comp in filter_pfas:
        PSDM_obj.test_range = np.array([PSDM_obj.k_data[comp]['ds']])
        comp, ds, ssqs, md, base = PSDM_obj.run_psdm_dsfit(comp)

        md[md<0.] = 0.
        md[md>data_df[comp][0]] = data_df[comp][0]
        if plot:
            plt.plot(md.index.values, md.values, label=comp)
        out_f = interp1d(md.index.values,\
                         md.transpose().values,\
                         fill_value='extrapolate')
        data_df[comp] = out_f(time_idx)[0]
        
    if plot:
        plt.legend(loc='center right')
        plt.ylabel('Concentration (ng/L)')
        plt.xlabel('Time (days)')
        plt.savefig('full_scale_'+PSDM_obj.carbon+'.png',dpi=300)
        plt.close()
    
    data_df[data_df<0]=0. #resets negatives to zero
    
    writer = pd.ExcelWriter(PSDM_obj.project_name+'_'+PSDM_obj.carbon+'.xlsx') 
    data_df.to_excel(writer, 'model_fit')
    writer.save()

    small_rotation = total_beds%beds_per_cycle
    
    if small_rotation == 0:
        #all cycles the same
        weights = int(total_beds/beds_per_cycle) * [float(beds_per_cycle)/total_beds]
    else:
        weights = [float(small_rotation)/total_beds] + \
                  int(total_beds/beds_per_cycle)*[float(beds_per_cycle)/total_beds]
        # having small rotation be first, is the worst case scenario
        # it means the smallest percentage of beds is new
        # having small rotation last is the best case scenario, not used.
    
    summed = data_df.transpose().sum()
    
    min_cycle = 1 #days
    
    num_cycles = np.ceil(float(total_beds)/beds_per_cycle)
    if num_cycles*beds_per_cycle > total_beds:
        print('number of beds per cycle may result in non-uniform cycling.\
              assuming new beds have fewest numbers')
    
    bed_info = []
    count = 0
    bed_c = 1
    for bed in range(total_beds):
        if count <= beds_per_cycle:
           bed_info.append(bed_c)
        count+=1
        if count == beds_per_cycle:
            count = 0
            bed_c+=1
           
    function = interp1d(summed.index,summed.values,fill_value='extrapolate')
    aa = np.arange(5*PSDM_obj.max_days)
    summed = pd.Series(function(aa),index = aa)
    
    best_cycle = np.zeros(min_cycle)
    best_val = min_cycle*1
    try:
        for i in range(min_cycle, PSDM_obj.max_days,1):
            tmp = np.zeros(i)
            for j in range(max(bed_info)):
                tmp += weights[j]*(summed[(summed.index>=(j*i))&(summed.index<((j+1)*i))].values)
            if tmp.max() <= target_conc:
                best_cycle = tmp*1.
                best_val = i
            else:
                break
    except Exception:
        best_val = PSDM_obj.max_days
        best_cycle = np.zeros(PSDM_obj.max_days)

    #return values to original values
    PSDM_obj.optimize_flag = opt_flg
    PSDM_obj.compounds = compounds
    PSDM_obj.test_range = test_range
    PSDM_obj.xdata = xdata
    
    return best_val, best_cycle


# =============================================================================
#     ANALYSIS FEATURES
# =============================================================================
def analyze_all(PSDM_obj):
    pool = mp.Pool(processes = PSDM_obj.processes)
    runs = [[i] for i in PSDM_obj.compounds]
    results = pool.starmap_async(PSDM_obj.run_psdm_kfit, runs)
    #runs all available compounds
    pool.close()
    
    real_results = results.get()
    PSDM_obj.real_results = real_results * 1
    for i in real_results:
        comp, k, xn, ssqs, md = i
        PSDM_obj.k_data[comp]['K'] = k
        PSDM_obj.k_data[comp]['1/n'] = xn
        
        #makes the maximum of the plot the 25% percentile value, to better
        #highlight minimum wells
        if PSDM_obj.optimize_flag:
            ssqs[ssqs>=np.percentile(ssqs.values,15)] = np.percentile(ssqs.values, 15)
            
            ##### plot the ssq space
            plt.figure()
            plt.contourf(ssqs.columns.values, ssqs.index.values,\
                         ssqs.values)
            
            min_val = find_minimum_df(ssqs)
  
            best_val_xn = min_val.columns
            best_val_k = min_val.index
            plt.plot(best_val_xn, best_val_k, 'rx')
            plt.title(comp+' - '+PSDM_obj.carbon)
            plt.xlabel('1/n')
            plt.ylabel('K-multiplier')
            plt.savefig(comp+'-'+PSDM_obj.carbon+'.png',dpi=300)
            plt.close()
        
        #plot the best fit
        dates = PSDM_obj.data_df.index
    
        plt.plot(dates, PSDM_obj.data_df[PSDM_obj.influent][comp], 
                 marker='.', ls=':', color='gray', label='Influent')
        plt.plot(dates, PSDM_obj.data_df[PSDM_obj.carbon][comp], 
                 marker='+', ls='None', color='black', label='Effluent')
        plt.plot(md.index, md.values, color='black', label='PSDM')
        plt.legend()
        plt.title(comp+' - '+PSDM_obj.carbon+'\n'+\
                  'K='+repr(round(k,3))+'   1/n='+\
                  repr(round(xn,3)))
        plt.xlabel('Time (days)')
        plt.ylabel('Concentration (ng/L)')
        plt.savefig(comp+'-'+PSDM_obj.carbon+'_model.png',dpi=300)
        plt.close()
        
        ti.sleep(.1)
        plt.close('all')

def analyze_all_ds(PSDM_obj):
    pool = mp.Pool(processes = PSDM_obj.processes)
    runs = [[i] for i in PSDM_obj.compounds]
    results = pool.starmap_async(PSDM_obj.run_psdm_dsfit, runs)
    #runs all available compounds
    pool.close()
    
    #create a new row for ds data
    PSDM_obj.k_data.loc['ds'] = PSDM_obj.k_data.loc['K'] * 1.
    PSDM_obj.k_data.loc['base_ds'] = PSDM_obj.k_data.loc['K'] * 1.
    
    real_results = results.get()
    for i in real_results:
        comp, ds, ssqs, md, base = i
        PSDM_obj.k_data[comp]['ds'] = ds * 1.
        PSDM_obj.k_data[comp]['base_ds'] = base * 1.
        
        #makes the maximum of the plot the 25% percentile value, to better
        #highlight minimum wells
        if PSDM_obj.optimize_flag:
            
            ##### plot the ssq space
            plt.figure()
            plt.plot(ssqs.index, ssqs.values)
            
            min_val = ssqs[ssqs==ssqs.min()]
  
            best_val_ds = min_val.index[0]
            plt.title(comp+' - '+PSDM_obj.carbon)
            plt.xlabel('Ds')
            plt.xscale('log')
            plt.ylabel('ssq')
            plt.savefig(comp+'-'+PSDM_obj.carbon+'.png',dpi=300)
            plt.close()
        
        #plot the best fit
        dates = PSDM_obj.data_df.index
    
        plt.plot(dates, PSDM_obj.data_df[PSDM_obj.influent][comp], marker='x',\
                 label='Influent')
        plt.plot(dates, PSDM_obj.data_df[PSDM_obj.carbon][comp], marker='o',\
                 label='Effluent')
        plt.plot(md.index, md.values, label='PSDM')
        plt.legend()
        plt.title(comp+' - '+PSDM_obj.carbon+'\n'+\
                  'K='+repr(round(PSDM_obj.k_data[comp]['K'],3))+'   1/n='+\
                  repr(round(PSDM_obj.k_data[comp]['1/n'],3))+'\nDs='+\
                  '{:0.4e}'.format(best_val_ds * base))
        plt.xlabel('Time (days)')
        plt.ylabel('Concentration (ng/L)')
        plt.savefig(comp+'-'+PSDM_obj.carbon+'_model.png',dpi=300)
        plt.close()
        
        ti.sleep(.1)
        plt.close('all')