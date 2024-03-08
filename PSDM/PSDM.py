# -*- coding: utf-8 -*-
"""
Combined GAC PSDM Model

Recoded from FORTRAN - AdDesignS - PSDM Model
Developed by:
David R. Hokanson, David W. Hand, John Crittenden, Tony N. Rogers, Eric J. Oman

Python Version
@author: Jonathan Burkhardt @UCChEJBB
         Levi Haupert
         
         
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
import warnings
warnings.simplefilter("ignore")

from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# import scipy as sp
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import interp1d
from scipy.stats import linregress
from scipy.optimize import curve_fit
import multiprocessing as mp
import time as ti #time as a variable is used in code, so ti is used

#Read in all associated PSDM functions
from PSDM_functions import min_per_day, lpg, spar_Jac, foul_params, kf_calc
from PSDM_functions import find_minimum_df, tortuosity, calc_solver_matrix
from PSDM_functions import density, viscosity, recalc_k, generate_grid
from PSDM_functions import interp 
from PSDM_functions import process_input_data, process_input_file
from PSDM_functions import logistic, filter_compounds

def run_MP_helper(test_column, k, invN, compound, k_mult):
    mp.freeze_support()
    
    test_column.test_range = np.array([k])
    test_column.xn_range = np.array([invN])
    try: # if it fails, just make a huge ssq. 
        compound, k, xn, ssqs, model_data = test_column.run_psdm_kfit(compound)
    except:
        ssqs = pd.DataFrame(1e9, columns=[invN], index=[k_mult] )
    return [k, invN, ssqs.values[0][0], compound, k_mult]

class PSDM():
    def __init__(self, column_data, comp_data, rawdata_df, **kw):
        '''
        column_data: contains specific characteristics of the column (pandas_df)
            'L'     = length (cm)
            'diam'  = diameter (cm)
            'wt'    = mass of carbon (g)
            'flrt'  = flow rate (units/min) see flow_type for units
            'rhop'  = apparent density (gm/ml)
            'rad'   = particle radius (cm) 
            'tortu' = tortuosity
            'psdfr' = pore to surface diffusion ratio
            'influentID' = influent identifier (text)
            'effluentID' = effluent identifier (text)
            'units' = may have units also specified
            --provided/processed by 'process_input_file' function
            
        comp_data: dataframe with compound information
            MW, MolarVol, BP, Density, Solubility (unused), VaporPress (unused)
        data_df: contains raw data associated with the column
            should contain influent and effluent data for the column
        
        keywords:
            project_name: (default = PSDM)
            nr: number of radial collocation points (int, default = 14)
            nz: number of axial collocation points (int, default = 19)
            1/n: Freundlich 1/n value (default = 0.45)
            
            temp: temperature in Celsius (default = 20)
            
            time_type: hours, days, mins (default = 'days')
            flow_type: gal, ml, L (default = 'ml')
            conc_type: ug, ng (per liter) (default = 'ng')
            
            break_type: calc, force (default = 'calc')
            brk_df: breakthrough dataframe (default = None)
            k_data: preprocessed k_data DataFrame (default = None)
            
            duration: length of simulation (default = max of data provided)
            
            water_type:
                default = 'Organic Free'
            chem_type:
                default = 'halogenated alkenes'
            
            mass_transfer:
                can be a dictionary with "{compound: {'kf': value, 'dp': value, 'ds':value}}"
                or
                can be a pandas dataframe with columns of compounds and index of ['kf','dp', 'ds']
                
                A value of zero will indicate that a user wants the correlation value.
                for dictionary, all values of 'kf', 'dp', or 'ds' do not have to be provided.
            
            test_range: default = np.linspace(1, 5, 41)
            xn_range:   default = np.arange(0.20, 0.95, 0.05)
            
            
        '''
        self.project_name = kw.get('project_name','PSDM')
        
        #collocation initilazation
        self.nc = kw.get('nr', 14) #set number of radial points, or 8
        self.mc = kw.get('nz', 19) #set number of axial points, or 11
        self.nz = self.mc * 1
        self.ne = kw.get('ne', 1) # number of finite elements
        solver_data = calc_solver_matrix(self.nc, self.mc, self.ne)
        self.wr = solver_data['wr']
        self.az = solver_data['az']
        self.br = solver_data['br']
        if self.mc != solver_data['mc']:
            ''' corrects for OCFE change to nz'''
            self.mc = solver_data['mc']
        self.nd = self.nc - 1
       
        #set up temperature dependant values
        self.temp = kw.get('temp', 20)
        self.vw = viscosity(self.temp)
        self.dw = density(self.temp)
        
        #define unit conversions (time)
        if 'time' not in column_data.index:
            self.time_type = kw.get('time_type', 'days') 
            if self.time_type == 'days': #base units in minutes
                self.t_mult = min_per_day
            elif self.time_type == 'hours':
                self.t_mult = 60. ## min_per_hour
            else:
                self.t_mult = 1. ## assumes minutes
        else:
            self.time_type = column_data.loc['time']
            self.t_mult = column_data.loc['t_mult']
        
        # flow units    
        if 'flow_type' not in column_data.index:
            self.flow_type = kw.get('flow_type', 'ml') 
            #deal with flow, conv. to liters
            if self.flow_type == 'gal':
                self.flow_mult = lpg
            elif self.flow_type == 'ml':
                self.flow_mult = 1e-3
            else:
                self.flow_mult = 1.   ## assumes liters    
        else:
            self.flow_type = column_data.loc['flow_type']
            self.flow_mult = column_data.loc['flow_mult']
        
        # concentration units
        if 'units' not in column_data.index:
            self.conc_type = kw.get('conc_type', 'ng') 
            #deal with mass, conv. to ug
            if self.conc_type == 'ug':
                self.mass_mul = 1.
            elif self.conc_type == 'ng':
                self.mass_mul = 1e-3
            elif self.conc_type == 'mg':
                self.mass_mul = 1e3
            else:
                print('conc_type is invalid, mg/ug/ng are valid options')
                self.mass_mul = 1.        ## assumes ug provided
        else:
            self.conc_type = column_data.loc['units']
            self.mass_mul = column_data.loc['mass_mul']
        
        #initialize column characteristics #process dataframe
        self.L = column_data['L']
        self.diam = column_data['diam']
        self.wt = column_data['wt']
        self.flrt = column_data['flrt']# * self.flow_mult
        self.rhop = column_data['rhop']
        self.rhof = column_data['rhof']
        self.rad = column_data['rad']
        self.tortu = column_data['tortu']
        self.psdfr = column_data['psdfr']
        self.epor = column_data['epor']
        self.influent = column_data['influentID']
        self.carbon = column_data.name
        self.duration = kw.get('duration',\
                               np.max(rawdata_df.index.values))
        #calculate other fixed values
        self.area = np.pi*(self.diam**2)/4.
        self.bedvol = self.area * self.L
        self.ebed = 1. - self.wt/(self.bedvol*self.rhop)
        self.tau = self.bedvol * self.ebed * 60./self.flrt
        self.sf = 0.245423867471 * self.flrt/self.area # gpm/ft**2
        self.vs = self.flrt/(60.*self.area)
        self.re = (2.*self.rad*self.vs*self.dw)/(self.ebed*self.vw)
        
        #calculate Empty Bed Contact Time (EBCT)
        self.ebct = self.area*self.L/self.flrt 
        
        #data
        self.data_df = rawdata_df
        self.comp_df = comp_data #information about the compounds
        self.xn = kw.get('xn', 0.45)
        self.xdata = rawdata_df.index.values
        
        #breakthrough type, data
        self.brk_type = kw.get('break_type','calc')
        self.brk_df = kw.get('brk_df', None)
        if self.brk_type == 'force': # self.brk_df != None and 
            self.brk_df['breakday'][self.brk_df['breakday']=='assume']=3000
            self.brk_df['breakday'][self.brk_df['breakday']=='est']=2000
            self.brk_df['breakday'][self.brk_df['breakday']=='no']=1000
            for i in self.brk_df.index:
                self.brk_df.iloc[i]['breakday'] = \
                int(self.brk_df.iloc[i]['breakday']) #convert everything to interger
                
        #filter
        self.compounds = kw.get('compounds',
                                list(rawdata_df.columns.levels[1]))
        self.num_comps = len(self.compounds) # used for multi_comp
        self.__y0shape = (self.num_comps, self.nc+1, self.mc)
        self.__altshape = (self.num_comps, self.nc, self.mc)
        
        # precalculate jacobian sparsity matrix
        self.jac_sparse = spar_Jac(self.num_comps, self.nc, self.nz, self.ne)
        self.jac_sparse_single = spar_Jac(1, self.nc, self.nz, self.ne)
             
        self.water_type = kw.get('water_type','Organic Free')
        self.chem_type = kw.get('chem_type', 'halogenated alkenes')
        
        self.test_range = kw.get('test_range', np.linspace(1, 5, 41))
        self.xn_range = kw.get('xn_range', np.arange(0.20, 0.95, 0.05))
        
        #handling for multiprocessing
        self.processes = kw.get('mp', mp.cpu_count())
        
        self.optimize_flag = kw.get('optimize', False)
        if len(self.test_range)==1 and len(self.xn_range)==1:
            self.optimize_flag = False
        
        #set the fouling parameters
        #expected that rk1-rk4 are pandas.Series
        self.rk1_store, self.rk2_store, self.rk3_store, self.rk4_store = self.__get_fouling_params()
        
        #set maximum cycles to consider for full scale
        self.max_days = kw.get('max_days', 2000) #number of days, 2000 default
        self.timeout = kw.get('timeout', 300) # 30 second default timeout
        self.plot_output = kw.get('plot', False)
        self.file_output = kw.get('file_out', True) #need to turn on
        self.solver = kw.get('solver', 'BDF') #defaults to BDF
        
        self.fouling_dict = {} #set up empty dictionary for fouling functions
        
        max_time = np.max([self.max_days, self.duration]) #doesn't do anything at the moment for k_fit
        
        self.time_vals = np.linspace(0, max_time*self.t_mult, 500)
        self.fouling_dict= self.__fouled_k_new(self.time_vals) ## removed self.k_data.loc['K'], 
        
        # calculate initial values 
        # might need to add more here
        k_data = kw.get('k_data', [])
        self.k_data_input_type = type(k_data)

        ## used in estimating K used in functions
        self.k_by_xn_factor = {} ## relationship between K and 1/n for a calculated q, should be dictionary of interpolating functions
        self.foul_mult_estimates = {comp: 1 for comp in self.compounds}  # creates list with 1x factor stored ## used to store estimate of fouling impact
        
        if len(k_data) == 0:
            self.k_data = pd.DataFrame(index=['K','1/n', 'q', 'brk','AveC'], \
                                        columns=self.compounds)
            for comp in self.compounds:
                k, q, classifier, brk, aveC, k_func, foul_mult_est = self.__calculate_capacity(comp)
                self.k_data[comp] = np.array([k, self.xn, q, brk, aveC])
                
                self.k_by_xn_factor[comp] = k_func
                self.foul_mult_estimates[comp] = foul_mult_est
                
        else:
            self.k_data = k_data
            if self.brk_type=='force':
                tmp_assume = self.brk_df[(self.brk_df['carbon']==self.carbon)&\
                                         (self.brk_df['breakday']!=3000)]
                impacted = tmp_assume['compound'].values
                for comp in self.compounds:
                    if comp in impacted:
                        k, q, classifier, brk, aveC, k_func, foul_mult_est = self.__calculate_capacity(comp)
                        self.k_data[comp] = np.array([k, self.xn, q, brk, aveC])
                        
                        self.k_by_xn_factor[comp] = k_func
                        
                        self.foul_mult_estimates[comp] = foul_mult_est
        
        self.k_data_orig = self.k_data
        
        if False:
            plt.figure()
            for key in self.fouling_dict.keys():
                plt.plot(self.time_vals/self.t_mult/7, \
                         self.fouling_dict[key](self.time_vals),\
                         label=key)
            plt.legend()
        

        ### adding mass transfer keyword
        self.mass_transfer = pd.DataFrame(0.0, columns=self.compounds, index=['kf','dp','ds'])
        mass_transfer_data = kw.get('mass_transfer', 'None')
        if type(mass_transfer_data) == str:
            pass
        elif type(mass_transfer_data) == pd.core.frame.DataFrame:
            try:
                for comp in mass_transfer_data.columns:
                    for val_type in mass_transfer_data.index:
                        self.mass_transfer.loc[val_type.lower(), comp] = mass_transfer_data[comp][val_type]
            except:
                print('Error in inputing mass_transfer values')
        elif type(mass_transfer_data) == dict:
            try:
                for comp in mass_transfer_data.keys():
                    for val_type in mass_transfer_data[comp].keys():
                        self.mass_transfer.loc[val_type.lower(), comp] = mass_transfer_data[comp][val_type]
            except:
                print('Error in inputing mass_transfer values')
        
        else:
            print('mass_transfer input must be a dictionary: {compound: {kf|dp|ds: value}}')
        
        self.__set_backups() ### set up backups

# =============================================================================
# end __init__
# =============================================================================
    def __get_fouling_params(self):
        '''
        water= [Rhine, Portage, Karlsruhe, Wausau, Haughton]
        chemical=  [halogenated alkanes, halogenated alkenes, trihalo-methanes
                    aromatics, nitro compounds, chlorinated compounds, phenols
                    PNAs, pesticides]
        '''
        a1, a2, a3, a4 = foul_params['water'][self.water_type]
        if self.chem_type != 'PFAS':
            b1, b2 = foul_params['chemical'][self.chem_type]
            dummy_array = np.ones(len(self.compounds))
            b1 = pd.Series(b1 * dummy_array, index=self.compounds)
            b2 = pd.Series(b2 * dummy_array, index=self.compounds)
        else:
            pfas_dict = foul_params['chemical']['PFAS']
            dict_keys = pfas_dict.keys() #available [a, b] pairs
            b1 = pd.Series(index=self.compounds) #start empty storage pd.Series
            b2 = pd.Series(index=self.compounds)
            
            for comp in self.compounds:
                if comp in dict_keys:
                    b1[comp] = pfas_dict[comp][0]
                    b2[comp] = pfas_dict[comp][1]
                    #forces the intercept of molar K reduction to 100% at t=0
                else: #if no specific value provided, return 'Average'
                    b1[comp] = pfas_dict['Ave'][0]
                    b2[comp] = pfas_dict['Ave'][1]
                    
        rk1 = b1 * a1 + b2
        rk2 = b1 * a2
        rk3 = b1 * a3
        rk4 = b1 * a4 #no factor of 100, in exponent (should there be b1?)
        return rk1, rk2, rk3, rk4
    
    def __fouled_k_new(self, t):
        # works on unconverted time
        # only works on multiplier
        if type(t) == np.ndarray:
            data_store = {}
            for comp in self.compounds:
                k_mult_pd = self.rk1_store[comp] + \
                            self.rk2_store[comp] * t + \
                            self.rk3_store[comp] * np.exp(self.rk4_store[comp] * t)
                
                k_mult_pd[k_mult_pd < 1e-3] = 1e-3
                
                data_store[comp] = interp1d(t,\
                                            k_mult_pd, \
                                            fill_value='extrapolate')
            return data_store
        
    def __calculate_capacity(self, compound):
        ## retunrs k, q_meas, breakthrough_code, breakthrough_time, aveC, k_function, foul_mult_est
        k = 0.
        q_meas = 0.
        breakthrough_code = 'none'
        breakthrough_time = self.duration
        foul_mult_est = 1 ### estimates multiplier related to fouling
        
        flow_per_day = self.flrt * self.flow_mult * self.t_mult ## should return L/day
        carbon_mass = self.wt   ## should be in grams
                
        xn_f_range = np.arange(0.15, 1.001, 0.01)  # sets up xn_range for returned interpolating function
        times_to_test = np.arange(self.duration + 1) ## returns a list of days 
        
        k_function = interp1d(xn_f_range, np.ones(len(xn_f_range))) ## creates empty function so something will be returned
        
        ### get influent and effluent data
        infl = self.data_df[self.influent, compound]
        infl[infl == 0] = 1e-3  ### prevents divide by zero error
        aveC = infl.mean()              #calculates average influent concentration
        
        effl = self.data_df[self.carbon, compound]
        
        ## create interpolating functions for influent and effluent
        f_inf = interp1d(infl.index, infl.values, fill_value='extrapolate')
        f_eff = interp1d(effl.index, effl.values, fill_value='extrapolate')
        
        xdata = self.xdata   ## can maybe get rid of
        
        donothing = False #set trigger for later process to False
        if infl.sum() == 0:
            print(f'No mass in influent for {compound}. Skipping.')
            donothing = True
              
        brk_found = False
        if self.brk_type == 'calc' and not donothing:
            if effl.sum() == 0:  ## no Effluent breakthrough at all
                print(f'Insufficient data to estimate breakthrough for {compound}. Returning minimum K estimate.')
                breakthrough_code = 'minimum'
                brk_found = True
                
                ## likely remove. duplicative
                # infl_load, _ = quad(f_inf, 0, breakthrough_time)
                # q = (infl_load) * flow_per_day * self.mass_mul / carbon_mass
                
                # aveC = np.mean(f_inf(times_to_test[times_to_test <= breakthrough_time]))
                
                # k = q / (aveC * self.mass_mul) ** self.xn
                
                # k_s = q / (aveC * self.mass_mul) ** xn_f_range
                # k_function = interp1d(xn_f_range, k_s, fill_value='extrapolate')
            elif np.count_nonzero(effl.values) == 1:
                ## if only one data point exceeds 0, may be able to estimate 
                if effl.iloc[-1].values[0] > 0 and effl.iloc[-1].values[0]/aveC >= 0.25:
                    ### if last value exceeds 25% of aveC, estimate line from last and second to last points???
                    ### hope logistic function finds solution?
                    pass

                else:
                    ## last point is zero or less than 25% of average C, means non-zero value is a middle value or unreliable way to estimate from this
                    brk_found = True
                    print(f'Insufficient data to estimate breakthrough for {compound}. Returning minimum K estimate.')
                    breakthrough_code = 'minimum'
                    
            if effl.max() >= infl.min() and not brk_found:  ## suggests likely crossover/breakthrough
                
                for t_idx in times_to_test: ## searches days to find cross_over
                    if not brk_found:
                        if f_eff(t_idx) >= f_inf(t_idx):
                            
                            breakthrough_time = (t_idx - 1) * 1 #shifts breakthrough to previous day
                            brk_found = True
                           
                            breakthrough_code = 'breakthrough'
                
                if not brk_found:
                    ### if previous search fails, check if C > aveC
                    for t_idx in times_to_test:
                        if not brk_found:
                            if f_eff(t_idx) >= aveC: ### must exceed 5% above average C. Replaces correlation finding.
                                ### may cause issues with highly variable influent
                                breakthrough_time = (t_idx - 1) * 1 #shifts breakthrough to previous day
                                brk_found = True

                                breakthrough_code = 'exceed_aveC, implied'
                    
            if not brk_found:
                ## Try fitting to a simplified logistic function first
                c_bound_low = np.min([0.99 * aveC, effl.values[0]])
                c_bound_high = np.max([aveC, effl.values[0]])
                try:
                    params, pcov = curve_fit(logistic, effl.index, effl.values, 
                                             p0=(aveC, 5, 0.1),
                                             bounds=((c_bound_low, 0, 0),
                                                     (c_bound_high, 1e3, 1)),
                                             maxfev=10000
                                             )
                    
                    breakthrough_time = np.round((params[1] + np.log(99)) / params[2], 0) ## rounds to nearest day, for 99% breakthrough relative to C0/AveC
                    
                    if breakthrough_time >= self.duration: ## rejects solution if breakthrough time doesn't exceed data duration (should have been caught)
                        brk_found = True
                        
                        breakthrough_code = 'logisitc'
                        
                        ## add fictitous point at end for integration step
                        infl.loc[breakthrough_time] = aveC
                        effl.loc[breakthrough_time] = logistic(breakthrough_time, *params)
                        
                        # update interpolating functions
                        f_inf = interp1d(infl.index, infl.values, fill_value='extrapolate')
                        f_eff = interp1d(effl.index, effl.values, fill_value='extrapolate')
                    elif breakthrough_time > 0.25 * self.duration:
                        ### may need to change, checks if predicted breakthrough point not at very beginning
                        brk_found = True
                        breakthrough_code = 'logistic'

                        
                except Exception as e:
                    print('Logistic search failed: ', e)
                    
                if not brk_found:
                    ### Try linearization as final check
                    try:
                        ### Try linearization
                        if effl.values[-1] > 0 and effl.values[-1]/aveC >= 0.25:
                            ### if last value exceeds 25% of aveC, estimate line from last and second to last points???
                            breakthrough_code = 'linear, single point'
                            brk_found = True
                            
                            slope = effl.values[-1]/(effl.index[-1] - effl.index[-2])
                            intercept = - effl.values[-1] * effl.index[-2] / (effl.index[-1] - effl.index[-2])
                            
                        else:
                            ### use linearization function
                            yte = effl.values  ## convenience variable
                            
                            nZc = np.count_nonzero(yte)
                            
                            if nZc < len(yte) and yte[-(nZc+1)] != 0: ## makes sure there is enough data to find non-zero before
                                nZc += 1 ### adds another point to the linearization search area
                            x = xdata[-(nZc+1):].astype('float64')
                            y = yte[-(nZc+1):].astype('float64')
                            
                            slope, intercept, *extra = linregress(x, y)
                            breakthrough_code = 'linear regression'
                            brk_found = True
                            
                        breakthrough_time = np.round((aveC - intercept)/slope,0)
                        
                        ## add fictitous point at end for integration step
                        infl.loc[breakthrough_time] = aveC
                        effl.loc[breakthrough_time] = aveC
                        
                        # update interpolating functions
                        f_inf = interp1d(infl.index, infl.values, fill_value='extrapolate')
                        f_eff = interp1d(effl.index, effl.values, fill_value='extrapolate')
                    
                    except Exception as e:
                        print('Linearization failed: ', e)
                        
            
            if brk_found:  ## process breakthrough estimate
                infl_load, _ = quad(f_inf, 0, breakthrough_time) ## influent loading, ignores error
                effl_rem, _ = quad(f_eff, 0, breakthrough_time) ## effluent removal, ignores error
        
                ## calculate q, ug/g
                q_meas = (infl_load - effl_rem) * flow_per_day * self.mass_mul / carbon_mass
                
                aveC = np.mean(f_inf(times_to_test[times_to_test <= breakthrough_time]))
                
                k = q_meas / (aveC * self.mass_mul) ** self.xn
                
                k_s = q_meas / (aveC * self.mass_mul) ** xn_f_range
                k_function = interp1d(xn_f_range, k_s, fill_value='extrapolate')
            
            else:
                print(f'WARNING: No Capacity Estimated for {compound}. Assuming data duration for capacity.')

                infl_load, _ = quad(f_inf, 0, breakthrough_time) ## influent loading, ignores error
                effl_rem, _ = quad(f_eff, 0, breakthrough_time) ## effluent removal, ignores error
        
                ## calculate q, ug/g
                q_meas = (infl_load - effl_rem) * flow_per_day * self.mass_mul / carbon_mass
                
                aveC = np.mean(f_inf(times_to_test[times_to_test <= breakthrough_time]))
                
                k = q_meas / (aveC * self.mass_mul) ** self.xn
                
                k_s = q_meas / (aveC * self.mass_mul) ** xn_f_range
                k_function = interp1d(xn_f_range, k_s, fill_value='extrapolate')


            
                
        elif self.brk_type == 'force':# and self.brk_df != None:
            brk_df = self.brk_df
            maxx = self.duration
            brkdy = brk_df[(brk_df['carbon']==self.carbon) & \
                           (brk_df['compound'] == compound)]['breakday'].values[0]
            
            if brkdy < 2000:
                f_infl = interp1d(xdata, infl, fill_value='extrapolate')
                f_effl = interp1d(xdata, effl, fill_value='extrapolate')
                
                xdata = xdata[xdata <= brkdy]
                
                if brkdy == 1000: #report min, but no breakthrough
                    qtmp = quad(f_infl, 0, maxx, points=xdata)[0] - \
                           quad(f_effl, 0, maxx, points=xdata)[0]
                    aveC = np.mean(infl)
                else:
                    if brkdy in xdata:
                        aveC = quad(f_infl, 0, brkdy)[0]/brkdy
                    else:
                        tmpC = f_infl(brkdy)
                        tmpxdata = xdata[xdata<=brkdy]
                        tmpxdata = np.append(tmpxdata, tmpC)
                        aveC = quad(f_infl,0, brkdy)[0]/brkdy
                        xdata = np.append(xdata, brkdy)
                    qtmp = quad(f_infl, 0, brkdy, points=xdata)[0] - \
                           quad(f_effl, 0, brkdy, points=xdata)[0]
                    breakthrough_time = brkdy
                
                int_infl = quad(f_infl, 0, breakthrough_time, points=xdata)[0]
                int_effl = quad(f_effl, 0, breakthrough_time, points=xdata)[0]
                qtmp = flow_per_day * (int_infl - int_effl) * self.mass_mul
                q_meas = qtmp/carbon_mass #ug/g
                k = q_meas / ((aveC*self.mass_mul) ** self.xn) 
            
            elif brkdy == 2000: #estimate
                nZc = np.count_nonzero(effl)
                
                x = xdata[-(nZc+1):]
                y = effl[-(nZc+1):]
                A = np.vstack([x,np.ones(len(x))]).T
                m,c = np.linalg.lstsq(A,y)[0]   #fits line  

                if m > 0:
                    intersection = (infl.mean() - c)/m
                else:
                    intersection = maxx #aveC
                    
                breakthrough_time = intersection
                
                infl = np.append(infl, aveC)
                effl = np.append(effl, aveC)
                xdata = np.append(xdata, intersection)
                
                f_infl = interp1d(xdata, infl, fill_value='extrapolate')
                f_effl = interp1d(xdata, effl, fill_value='extrapolate')
                aveC = quad(f_infl, 0, breakthrough_time)[0]/\
                       breakthrough_time
                
                qtmp = quad(f_infl, 0, intersection)[0] - \
                       quad(f_effl, 0, intersection)[0]
                        
                q_meas = flow_per_day * qtmp * self.mass_mul /\
                         (carbon_mass) # ug/g
                k = q_meas / ((aveC * self.mass_mul) ** self.xn) 
                
            breakthrough_code = 'supplied'
        
        
        ## should return the averaged impact related to K reduction caused by fouling
        foul_mult_est = 1/np.mean(self.fouling_dict[compound](np.arange(breakthrough_time)*self.t_mult))
        
        # returns capacity in (ug/g)(L/ug)**(1/n)
        return k, q_meas, breakthrough_code, breakthrough_time, aveC, k_function, foul_mult_est

    def __set_backups(self):
        ## store initial values so they can be changed: backups flagged with _bup
        self.k_data_bup = self.k_data.copy()
        self.data_bup = self.data_df.copy()
        self.compounds_bup = self.compounds
        self.num_comps_bup = self.num_comps * 1
        self.jac_sparse_bup = self.jac_sparse
        self.y0shape_bup = self.__y0shape #(self.num_comps, self.nc+1, self.mc)
        self.altshape_bup = self.__altshape #(self.num_comps, self.nc, self.mc)
        self.mass_transfer_bup = self.mass_transfer.copy()
        
            ### column values to reset later
        self.L_bup = self.L * 1
        self.diam_bup = self.diam * 1 
        self.wt_bup = self.wt * 1 
        self.flrt_bup = self.flrt * 1 
        self.rhop_bup = self.rhop * 1 
        self.rhof_bup = self.rhof * 1 
        self.rad_bup = self.rad * 1
        self.tortu_bup = self.tortu * 1 
        self.psdfr_bup = self.psdfr * 1
        self.epor_bup = self.epor * 1

    def __reset_column_values(self):
        ### just a deep copy???
        # =============================================================================
        #         Reset values in column object
        # =============================================================================
        self.L = self.L_bup * 1
        self.diam = self.diam_bup * 1 
        self.wt = self.wt_bup * 1
        self.flrt = self.flrt_bup * 1
        self.rhop = self.rhop_bup * 1
        self.rhof = self.rhof_bup * 1
        self.rad = self.rad_bup * 1
        self.tortu = self.tortu_bup * 1 
        self.psdfr = self.psdfr_bup * 1
        self.epor = self.epor_bup * 1
        
        #calculate other fixed values
        self.area = np.pi*(self.diam**2)/4.
        self.bedvol = self.area * self.L
        self.ebed = 1. - self.wt/(self.bedvol*self.rhop)
        self.tau = self.bedvol * self.ebed * 60./self.flrt
        self.sf = 0.245423867471 * self.flrt/self.area # gpm/ft**2
        self.vs = self.flrt/(60.*self.area)
        self.re = (2.*self.rad*self.vs*self.dw)/(self.ebed*self.vw)
        
        # #calculate Empty Bed Contact Time (EBCT)
        self.ebct = self.area*self.L/self.flrt 
        
        self.k_data = self.k_data_bup.copy()
        self.data_df = self.data_bup.copy()
        self.compounds = self.compounds_bup
        self.num_comps = self.num_comps_bup * 1
        self.__y0shape = self.y0shape_bup 
        self.__altshape = self.altshape_bup
        self.jac_sparse = self.jac_sparse_bup #spar_Jac(self.num_comps, self.nc, self.nz, self.ne)
        self.mass_transfer = self.mass_transfer_bup.copy()  
        
        ## returns nothing, just resets class variables to originally calculated states for interative loop in model_uncertainty()

# =============================================================================
# End INIT and helper functions
# =============================================================================
    
    def run_all(self, plot=False, save_file=True, optimize='staged', 
                init_grid=5, init_loop=3):
        '''
        Parameters
        ----------
        plot : BOOL, optional
            Are plots generated during this function. The default is False.
        save_file : BOOL, optional
            Are results files generated during this function. The default is True.
        optimize : string, optional
            'brute' or 'staged' are acceptable options. The default is 'staged'.
        init_grid : INT, optional
            Number of grid points for staged optimizer. The default is 5.
        init_loop : INT, optional
            Number of refinement loops for staged optimizer. The default is 3.

        Current Behavior
        ----------------
        All files will be overwritten. No file checking is currently performed.
        
        Returns
        -------
        None. Optimized k & 1/n values are saved to Class object, 
        and can be saved in script that calls this function.

        '''
        ## TODO: Check this still functions with recent updates, 1/2024
        #forces single run to handle optimizing in this funciton, not run_psdm_kfit
        opt_flg = self.optimize_flag
        orig_test_range = self.test_range * 1.
        orig_xn_range = self.xn_range * 1.
        self.optimize_flag = False
        
        xn_rng_tmp = np.zeros(2)
        k_rng_tmp = np.zeros(2)
        des_xn = 0.025 # search accuracy for xn
        des_k = 0.1    # search accuracy for k multiplier
        
        if optimize == 'brute':
            init_loop = 0
            k_range = self.test_range
            xn_range = self.xn_range

        for compound in self.compounds:
            print(compound, ' running')

            k_val = self.k_data[compound]['K']
            q_val = self.k_data[compound]['q']
            
            inf = self.data_df[self.influent][compound]
            eff = self.data_df[self.carbon][compound]
            
            xn_range1 = orig_xn_range * 1.
            k_range1 = orig_test_range * 1.
            
            grid_num_xn = init_grid * 1
            grid_num_k = init_grid * 1
            
            for loop in range(init_loop + 1):
                data_store = [] #reset
                k_mult = {}
                models = {}

                if optimize == 'staged':
                    if loop == init_loop:
                        # expected different search space for final loop
                        xn_rng_tmp[0] = np.floor(xn_rng_tmp[0] / des_xn) 
                        xn_rng_tmp[1] = np.ceil(xn_rng_tmp[1] / des_xn)
                        
                        grid_num_xn = int(xn_rng_tmp.ptp() + 1)
                        xn_range1 = xn_rng_tmp * des_xn
                        if xn_range1[1] > orig_xn_range[-1]:
                            #corrects for overrun of xn_range
                            #may need to do the same for underrun of xn_rng
                            xn_range1[1] = orig_xn_range[-1]
                            grid_num_xn -= 1
                        
                        k_rng_tmp[0] = np.floor(k_rng_tmp[0] / des_k)
                        k_rng_tmp[1] = np.ceil(k_rng_tmp[1] / des_k)
                        grid_num_k = int(k_rng_tmp.ptp() + 1)
                        k_range1 = k_rng_tmp * des_k

                    xn_range, k_range = generate_grid(grid_num_xn,
                                                      grid_num_k,
                                                      loop_num=loop,
                                                      xn_range=xn_range1, 
                                                      k_range=k_range1)
                for xns in xn_range:
                    k_mult[xns] = recalc_k(k_val, q_val, self.xn, xns)
            
                for k in k_range:
                    models[k] = {}
                    for xn in xn_range:
                        self.test_range = np.array([k*k_mult[xn]])
                        self.xn_range = np.array([xn])
                        
                        comp, k_v, xn_v, ssqs, md = self.run_psdm_kfit(compound)
                        data_store.append([k, xn_v, ssqs.values[0][0]])
                        
                        models[k][xn] = md
                        
                data_pd = pd.DataFrame(data_store, columns=['K','1/n','ssq'])
                ssqs = pd.pivot_table(data_pd,values=['ssq'], 
                                                 index=['K'],
                                                 columns=['1/n'],
                                                 aggfunc=np.max,
                                                 )['ssq']
                min_val = find_minimum_df(ssqs)
                
                xn_grid = xn_range[1]-xn_range[0]
                k_grid = k_range[1]-k_range[0]
                
                min_xn = min_val.columns[0]
                min_k = min_val.index[0]
                idx_xn = xn_range.tolist().index(min_xn)
                idx_k = k_range.tolist().index(min_k)
                
                max_idx_xn = len(xn_range) - 1
                max_idx_k  = len(k_range) - 1
                
                if idx_xn == 0: #close to min xn side
                    xn_rng_tmp[0] = xn_range[0]
                    xn_rng_tmp[1] = xn_range[1]
                elif idx_xn == max_idx_xn: # close to max xn side
                    xn_rng_tmp[0] = xn_range[-2]
                    xn_rng_tmp[1] = xn_range[-1]
                else: #middle of search space
                    xn_rng_tmp[0] = xn_range[idx_xn-1]
                    xn_rng_tmp[1] = xn_range[idx_xn+1]
                    
                if idx_k == 0: #close to min k side
                    k_rng_tmp[0] = k_range[0]
                    k_rng_tmp[1] = k_range[1]
                elif idx_k == max_idx_k: # close to max k side
                    k_rng_tmp[0] = k_range[-2]
                    k_rng_tmp[1] = k_range[-1]
                else: #middle of search space
                    k_rng_tmp[0] = k_range[idx_k-1]
                    k_rng_tmp[1] = k_range[idx_k+1]
                
                if xn_grid < des_xn:
                    #can reduce search space
                    grid_num_xn = 3
                    if idx_xn == max_idx_xn:
                        grid_num_xn = 2
                        xn_rng_tmp[0] = xn_rng_tmp[1] - des_xn 
                    elif idx_xn == 0:
                        grid_num_xn = 2   
                        xn_rng_tmp[1] = xn_rng_tmp[0] + des_xn
                
                if k_grid < des_k:
                    #can reduce seach space
                    grid_num_k = 3
                    if idx_k == max_idx_k:
                        grid_num_k = 2
                        k_rng_tmp[0] = k_rng_tmp[1] - des_k 
                    elif idx_k == 0:
                        grid_num_k = 2
                        k_rng_tmp[1] = k_rng_tmp[0] + des_k
                    
                xn_range1 = xn_rng_tmp * 1.
                k_range1 = k_rng_tmp * 1.
                
            min_val = find_minimum_df(ssqs)
            best_val_xn = min_val.columns[0]
            best_val_k = min_val.index[0] * k_mult[best_val_xn]
            md = models[min_val.index[0]][best_val_xn]
            min_val = min_val.values[0][0]

            if plot:
                plt.figure()
                plt.plot(inf.index, inf.values, marker='x', label='influent')
                plt.plot(eff.index, eff.values, marker='.',
                         markerfacecolor='None',label='effluent')
                #plot model results
                plt.plot(md.index, 
                         md.values,
                         label=f"K: {best_val_k:.2f}, 1/n: {best_val_xn:.3f}")

                
                plt.legend()
                plt.title(f"{compound} - {self.carbon}")
                plt.savefig(f"{self.carbon}_{compound}.png",dpi=300)
                plt.close()
                
                plt.figure()
                ssqs[ssqs>np.percentile(ssqs, 25)] = np.percentile(ssqs, 25)
                plt.contourf(ssqs.columns, ssqs.index, ssqs.values)
                plt.title(f"{compound} - {self.carbon}")
                plt.savefig(f"ssq_{self.carbon}_{compound}.png",dpi=300)
                plt.close()
            
            if save_file:
                with pd.ExcelWriter(f"ssq_{self.carbon}-{compound}.xlsx") as writer:
                    ssqs.to_excel(writer, 'Sheet1')
                
                with pd.ExcelWriter(self.project_name+'_'+self.carbon + '-' + \
                                        compound+'.xlsx') as writer:
            
                    md.to_excel(writer, 'model_fit')
                    
                    inf.to_excel(writer, 'influent')
                    eff.to_excel(writer, 'effluent')
                    
                    data_tmp = pd.Series([self.re, best_val_k, best_val_xn,\
                                 min_val, self.ebct, self.sf], \
                                 index=['Re','K','1/n','ssq','ebct','sf'])
                    data_tmp.to_excel(writer, 'parameters')
            
            self.k_data[compound]['K'] = best_val_k
            self.k_data[compound]['1/n'] = best_val_xn

        self.xn_range = orig_xn_range
        self.test_range = orig_test_range
        self.optimize_flag = opt_flg #resets to original value
# end run_all()

    def run_all_smart(self, plot=False, 
                      save_file=True, file_name='PSDM_', 
                      pm=10, num=11, des_xn=0.025, 
                      search_limit=50):
        '''
        Smart Optimizer for K & 1/n fitting.
        Starts with estimates for effective fouling, and attemps to find path 
        to lowest ssq by following decreasing trends. 

        Given limited data, the resulting K & 1/n values may be close to the 
        "true" K & 1/n values, but should provide reasonable predictive value
        for modeling performance. model_uncertainty function can be used
        to predict sensitivity of model results to these paramaters. 
        
        
        NOTE: No checking of files is done by this script. If files exist
        with the same filenames as those generated by the script, they will be
        overwritten.

        Parameters
        ----------
        plot : BOOL, optional
            Flag for specifying whether to generate plots. The default is False.
        save_file : BOOL, optional
            Flag to designate if xlsx files should be saved. The default is True.
        file_name: string, optional
            Base name of files generated by this script. The default is True.
        pm : number, optional
            plus/minus range of additional test range for K (0-50). The default is 0.
            20 = 20% so additional loop of 0.8-1.2x Kmultiplier will be tested
            after 1/n is determined.
        num: INT, optional
            number of additional grids to search if pm != 0. The default is 11.
            Should be odd. This is used by np.linspace(1-pm, 1+pm, num)
        des_xn: float, optional
            defines grid size for 1/n search. The default is 0.025.
        search_limit: int or float, optional
            defines an upper limit for how many searches can occur. Default is 50.

        Returns
        -------
        None.
        Best K & 1/n values are stored to self.k_data

        '''
        def get_ssq(k, k_factor, xn, compound, ssq_storage):
            self.k_data.loc['K', compound] = k * k_factor * 1 
            self.k_data.loc['1/n', compound] = xn * 1 
            
            if np.isnan(ssq_storage.loc[np.round(k_factor, 6), np.round(xn,3)]):
                ## only runs if no ssq data is already available
                ## prevents duplicative run_psdm_kfit calls
                try:
                    _, _, _, ssqs, _ = self.run_psdm_kfit(compound)

                    # print(ssqs)
                
                    return ssqs.values[0][0]
                except Exception as e:
                    ## if error, return big number, won't become best_ssq
                    print(e)
                    return 1e20
            else:
                return ssq_storage.loc[np.round(k_factor,6), np.round(xn,3)]
            
        def get_base_k(compound, xn):
            if compound in self.k_by_xn_factor.keys():
                ## use precalcualted value from __calculate_capacity
                return self.k_by_xn_factor[compound](xn)
            else:
                ## return K from k_data
                return self.k_data[compound]['K']
        
        def get_fouling_factor(compound):
            if compound in self.foul_mult_estimates.keys():
                ## use precalcualted value from __calculate_capacity
                return self.foul_mult_estimates[compound]
            else:
                ## calculate estimate
                brk_day = self.k_data[compound]['brk']
                comp_fouling = self.fouling_dict[compound]
                integral, _ = quad(comp_fouling, 0, brk_day) #ignore error from quad
                return brk_day/integral

        ## store some initial conditions
        opt_flg = self.optimize_flag
        orig_test_range = self.test_range * 1.
        orig_xn_range = np.round(self.xn_range * 1.,5)
        orig_k_data_file_type = self.k_data_input_type
        self.k_data_input_type = int
        self.optimize_flag = False
        
        file_name = file_name.replace(' ', '_') #removes spacing for underscore
        if file_name[-1] != '-' or file_name[-1] != '_':
            file_name += '_' #adds spacer underscore to filename
            
        for compound in self.compounds:  ## search all compounds
            print(f'Running - {compound}')
            k_val, q_val, brk_day = self.k_data[compound][['K', 'q', 'brk']]
            
            ## set starting conditions
            best_xn = 0.45 ## reasonable starting point
            best_k_factor = get_fouling_factor(compound)
            base_k = get_base_k(compound, best_xn)
            best_k = base_k * best_k_factor 

            ## initially test increasing xn
            count = 0
            cases_tried = 0
            scale_k = best_k_factor * pm / 100 / num
            max_k = best_k_factor * (1 + pm/100)
            min_k = best_k_factor * (1 - pm/100)
            
            #starter variables
            xn = best_xn * 1
            k_factor = get_fouling_factor(compound)
            
            ssq_storage = pd.DataFrame(index=np.round(np.arange(min_k, max_k+scale_k/2, scale_k),6),
                                       columns=np.round(np.arange(0.2, 1+des_xn/2, des_xn),3))
                        
            best_ssq = get_ssq(base_k,  best_k_factor, best_xn, compound, ssq_storage)
            ssq_storage.loc[np.round(best_k_factor,6), np.round(best_xn,3)] = best_ssq * 1
            
            ## set starting index locations
            xn_idx = np.where(ssq_storage.columns == np.round(best_xn,3))[0][0]
            k_idx = np.where(ssq_storage.index == np.round(best_k_factor,6))[0][0]
            xn_idx_max = len(ssq_storage.columns) - 1 
            k_idx_max = len(ssq_storage.index) - 1

            ### attempt to get best starting location, assume initial guess at k is correct
            xn_search_idx = np.round(np.linspace(0, xn_idx_max, num=5),0).astype(int) ## should return 0, three middle points, end for rough search.
            
            for xn_test in xn_search_idx:
                xn = ssq_storage.columns[xn_test]
                k_factor = ssq_storage.index[k_idx]
                test_base_k = get_base_k(compound, xn)

                test_ssq = get_ssq(test_base_k, k_factor, xn, compound, ssq_storage)

                ssq_storage.loc[np.round(k_factor,6), np.round(xn,3)] = test_ssq * 1

                if test_ssq < best_ssq:
                    best_xn = xn * 1
                    best_k = k_factor * test_base_k * 1
                    best_ssq = test_ssq * 1

                    xn_idx = xn_test * 1 ### resets xn_idx to best starting location of group.

            ### search from best starting location
            shift_xn = 0
            shift_k = 0
            while cases_tried <= 5: 
                ### tries to seach in 4 directions. Uses this information to search again

                ## move central test location
                xn_idx += shift_xn
                k_idx += shift_k

                ### get current location ssq
                xn = ssq_storage.columns[xn_idx]
                k_factor = ssq_storage.index[k_idx]
                center_base_k = get_base_k(compound, xn)
                
                base_ssq = get_ssq(center_base_k,  k_factor, xn, compound, ssq_storage)
                if np.isnan(ssq_storage.loc[np.round(k_factor,6), np.round(xn,3)]):
                    ssq_storage.loc[np.round(k_factor,6), np.round(xn,3)] = base_ssq * 1
                
                if base_ssq < best_ssq:
                    best_xn = xn * 1
                    best_k = k_factor * center_base_k * 1
                    best_ssq = base_ssq * 1

                ### try shifting to the right (xn)
                if xn_idx < xn_idx_max:
                    right_xn = ssq_storage.columns[xn_idx + 1]
                    right_base_k = get_base_k(compound, right_xn)

                    right_ssq = get_ssq(right_base_k,  k_factor, right_xn, compound, ssq_storage)
                    if np.isnan(ssq_storage.loc[np.round(k_factor,6), np.round(right_xn,3)]):
                        ssq_storage.loc[np.round(k_factor,6), np.round(right_xn,3)] = right_ssq * 1

                    if right_ssq < best_ssq:
                        best_xn = right_xn * 1
                        best_k = k_factor * right_base_k * 1
                        best_ssq = right_ssq * 1
                else:
                    right_ssq = 1e20 ### right_ssq not available, make huge

                ### try shifting to the left
                if xn_idx >= 1:
                    left_xn = ssq_storage.columns[xn_idx - 1]
                    left_base_k = get_base_k(compound, left_xn)

                    left_ssq = get_ssq(left_base_k,  k_factor, left_xn, compound, ssq_storage)
                    if np.isnan(ssq_storage.loc[np.round(k_factor,6), np.round(left_xn,3)]):
                        ssq_storage.loc[np.round(k_factor,6), np.round(left_xn,3)] = left_ssq * 1

                    if left_ssq < best_ssq:
                        best_xn = left_xn * 1
                        best_k = k_factor * left_base_k * 1
                        best_ssq = left_ssq * 1
                else:
                    left_ssq = 1e20 ### left_ssq not available, make huge

                ## Decide where to move next (xn)
                if right_ssq <= base_ssq and xn_idx + 1 <= xn_idx_max:
                    shift_xn = 1
                elif left_ssq <= base_ssq and xn_idx - 1 >= 0:
                    shift_xn = -1
                else:
                    ## can't move in xn direction
                    shift_xn = 0

                ### try shifting up (k)
                if k_idx < k_idx_max:
                    upper_k_factor = ssq_storage.index[k_idx + 1]

                    upper_ssq = get_ssq(center_base_k,  upper_k_factor, xn, compound, ssq_storage)
                    if np.isnan(ssq_storage.loc[np.round(upper_k_factor,6), np.round(xn,3)]):
                        ssq_storage.loc[np.round(upper_k_factor,6), np.round(xn,3)] = upper_ssq * 1

                    if upper_ssq < best_ssq:
                        best_xn = xn * 1
                        best_k = upper_k_factor * center_base_k * 1
                        best_ssq = upper_ssq * 1
                else:
                    upper_ssq = 1e20 ### upper_ssq not available, make huge

                ### try shifting up (k)
                if k_idx >= 1:
                    lower_k_factor = ssq_storage.index[k_idx - 1]

                    lower_ssq = get_ssq(center_base_k,  lower_k_factor, xn, compound, ssq_storage)
                    if np.isnan(ssq_storage.loc[np.round(lower_k_factor,6), np.round(xn,3)]):
                        ssq_storage.loc[np.round(lower_k_factor,6), np.round(xn,3)] = lower_ssq * 1

                    if lower_ssq < best_ssq:
                        best_xn = xn * 1
                        best_k = lower_k_factor * center_base_k * 1
                        best_ssq = lower_ssq * 1
                else:
                    lower_ssq = 1e20 ### lower_ssq not availabale, make huge

                ## Decide where to move next (K)
                if upper_ssq <= base_ssq and k_idx + 1 <= k_idx_max:
                    shift_k = 1
                elif lower_ssq <= base_ssq and k_idx - 1 >= 0:
                    shift_k = -1
                else:
                    ## can't move in K direction
                    shift_k = 0
                

                ### if ready to end loop, save results
                if base_ssq == best_ssq:
                    ### if the central location is the best, break the loop
                    
                    self.ssq_storage = ssq_storage 
                    
                    # print('Best', best_xn, best_k)
                    self.k_data.loc['1/n', compound] = best_xn * 1
                    self.k_data.loc['K', compound] = best_k * 1

                    self.k_data_bup = self.k_data.copy() ### reset backup as well
                    cases_tried = 10

                elif shift_k == 0 and shift_xn == 0:
                    ## predicting no movement, break loop, 
                    ## may only get triggered on edges(?) but just in case
                    self.ssq_storage = ssq_storage 

                    # print('Best alt', best_xn, best_k)
                    self.k_data.loc['1/n', compound] = best_xn * 1
                    self.k_data.loc['K', compound] = best_k * 1

                    self.k_data_bup = self.k_data.copy() ### reset backup as well
                    
                    cases_tried = 10



        if save_file:
            with pd.ExcelWriter('best_fits-'+self.project_name+'.xlsx') as writer:
                self.k_data.to_excel(writer, 'Sheet1')            
                        
        #reset original values
        # self.optimize_flag = opt_flg
        self.test_range = orig_test_range * 1.
        self.xn_range = orig_xn_range * 1. 
        self.k_data_input_type = orig_k_data_file_type
        self.__reset_column_values()
      

#END RUN_ALL_SMART
    


# =============================================================================
# #### Begin run_psdm()
# =============================================================================

    def run_psdm(self):
        '''
        new run multi
        
        Allows multi component competitive adsorption, with fouling
        
        when number of components exceeds 2, and fouling is used:
            use nz = 3, and ne = ## (20-30) to control axial number.
            this provides more even spacing of axial points, which helps 
            reduce stiff transitions in ODE solver.
        '''
        #pull information out of solver_data/self
        wr = self.wr
        nc = self.nc
        mc = self.mc
        az = self.az
        br = self.br
        t_mult = self.t_mult
        vw = self.vw        #viscosity
        dw = self.dw        #density
        epor = self.epor
        rhop = self.rhop
        ebed = self.ebed
        tau = self.tau
        rad = self.rad
        compound_list = self.compounds
                
        water_type = self.water_type
        
        k_v = self.k_data.loc['K']
        xn_v = self.k_data.loc['1/n']
        
        mw = self.comp_df.loc['MW']           #molecular weight
        mol_vol = self.comp_df.loc['MolarVol'] # Molar volume
        
        if self.time_type == 'days':
            dstep = 0.25 * min_per_day
        else:
            dstep = 1.       
            
        #set up bindings for nonlocal varaibles
        cinf = 1.
        cout_f = 1.
        tconv = 1.
        time_dim = 1.
        time_dim2 = 1.
        ttol = 1.
        tstep = 1.
        ds_v = 1. 
                
        inf = self.data_df[self.influent]
        
        #convert cbo to molar values
        cbo = (inf * self.mass_mul / mw)[compound_list]
        time = (inf.index * t_mult).values # time data from index
        cb0 = cbo.loc[0]
        cb0[cb0 == 0.] = 1. # if influent at time = 0 is 0, reset to 1 
        cin = cbo/cb0 # convert to relative concentration
        
        tortu = self.tortu                             # tortuosity
        psdfr = self.psdfr                             # pore to surface diffusion ratio
        nd = nc - 1
        
        difl = 13.26e-5/(((vw * 100.)**1.14)*(mol_vol**0.589)) #vb
        sc = vw / (dw * difl)       #schmidt number
        
        #set film and pore diffusion
        multi_p = difl/(2*rad) # multiplier used for kf calculation
        ## calculate everything first, replace as needed
        kf_v = kf_calc(multi_p, self.re, sc, ebed, corr='Chern and Chien')
                    
        dp_v = (difl/(tortu))       #*column_prop.loc['epor'] #porosity not used in AdDesignS appendix, removed to match
        
        self.mass_transfer_data = self.mass_transfer.copy()
        
        def run():
            nonlocal cinf
            nonlocal cout_f
            nonlocal tconv
            nonlocal time_dim
            nonlocal time_dim2
            nonlocal ttol
            nonlocal tstep
            nonlocal ds_v # for output

            aau = np.zeros(mc)
            #==============================================================================
            # #converts K to (umole/g)(L/umole)**(1/n), assumes units of (ug/g)(L/ug)**(1/n)
            #==============================================================================
            molar_k = (k_v / mw / ((1. / mw)**xn_v))[self.compounds]

            xni = 1./xn_v
            
            ds_v = epor*difl*cb0*psdfr/(1e3*rhop*molar_k*cb0**xn_v)
            
            for cdx in self.mass_transfer.columns:
                if self.mass_transfer[cdx]['kf'] > 0.:
                    kf_v[cdx] = self.mass_transfer[cdx]['kf']
                if self.mass_transfer[cdx]['dp'] > 0.:
                    dp_v[cdx] = self.mass_transfer[cdx]['dp']
                if self.mass_transfer[cdx]['ds'] > 0.:
                    ds_v[cdx] = self.mass_transfer[cdx]['ds']
                    ### should this be before or after the fouling? 
            
            if water_type != 'Organic Free':
                ds_v /= 1e10  #original code =1e-30, but this caused instability
            
            ## save input values for later use
            self.mass_transfer_data.loc['kf'] = kf_v
            self.mass_transfer_data.loc['dp'] = dp_v
            self.mass_transfer_data.loc['ds'] = ds_v
            
            d = (ds_v/dp_v)[self.compounds]
            
            qe = molar_k * cb0**xn_v
            qte = qe.sum()
                       
            dgs = (rhop * qe * (1.-ebed) * 1000.)/(ebed * cb0)
            dgp = epor * (1. - ebed)/(ebed) 
            dg = dgs + dgp
            dgt = dg.sum()
            dg1 = 1. + dgt
            dgI = 1.0/dg
            edd = dgt/dg #dgt changed from dg1
            
            ym = qe/qte
            
            eds = (ds_v*dgs*tau/(rad**2))[self.compounds]
            eds[eds < 1e-130] = 1e-130 #set lower limit of eds
            edp = (dp_v*dgp*tau/(rad**2))[self.compounds]
                        
            # from orthog(n)
            beds = (eds + d*edp) * edd
            bedp = edp * (1. - d) * edd
            
            #depends on kf
            st = (kf_v * (1. -ebed) * tau/(ebed*rad))[self.compounds]
            stdv = st * dgt     # dgt changed from dg1
            
            #convert to dimensionless
            tconv = 60./(tau*dg1)
            self.tconv = tconv
            tstep = dstep * tconv
            ttol = time[-1] * tconv
            time_dim = time * tconv
            
            #set up time based influent data
            time_temp = np.arange(0, time[-1]+self.t_mult*10, 1)

            cinfA = np.zeros((self.num_comps, len(time_temp)))
            foulFA = np.zeros((self.num_comps, len(time_temp)))
            facA = np.zeros((self.num_comps, len(time_temp)))
            
            #precalculate time-based toruosity    
            if water_type != 'Organic Free':
                tortu = tortuosity(time_temp)
            else:
                tortu = np.ones(len(time_temp))
                
            #convert pandas to arrays, test for speed? 2-3x speed up per diffun
            ThreeDSize = (self.num_comps, 1, 1)
            TwoDSize = (self.num_comps, 1)

            ym_v = ym[compound_list].values
            ym_vA = ym_v.reshape(ThreeDSize) 
            xni_v = xni[compound_list].values
            xni_vA = xni_v.reshape(ThreeDSize)
            xn_vval = xn_v[compound_list].values
            xn_vA = xn_vval.reshape(ThreeDSize)
            beds_v = beds[compound_list].values
            beds_mult = np.multiply(np.tile(br, ThreeDSize), 
                                    beds_v.reshape(ThreeDSize))
            bedp_v = bedp[compound_list].values
            bedp_mult = np.multiply(np.tile(br, ThreeDSize),
                                    bedp_v.reshape(ThreeDSize))
            wr_A = np.tile(wr,TwoDSize)
            az_A = np.tile(az, ThreeDSize)
            stdv_v = stdv[compound_list].values
            dgI_v = dgI[compound_list].values
            stdv_dgI = stdv_v.reshape(TwoDSize)*\
                        dgI_v.reshape(TwoDSize)
            dgI_v = dgI_v.reshape(TwoDSize)
            xn_ymA = np.divide(xn_vA, ym_vA)
            
            for i in range(self.num_comps):
                comp = compound_list[i]
                cinf = interp1d(time,\
                                cin[comp].values, \
                                fill_value = 'extrapolate') 
                
                cinfA[i] = cinf(time_temp)
                facA[i] = (1./tortu-d[comp])/(1.-d[comp])
                foulFA[i] = (1./self.fouling_dict[comp](time_temp))**xni[comp]
            
            cinfA = cinfA.transpose()
            foulFA = foulFA.transpose()
            facA = facA.transpose()
            
            #initialize storage arrays/matrices
            n = self.num_comps * (nc+1)*mc
            y0 = np.zeros(n)
            # yt0 = np.zeros(self.__altshape)
            # z = np.zeros(self.__altshape)
            # q0 = np.zeros(self.__altshape)
            aau = np.zeros((self.num_comps, mc))
            aau2 = np.zeros((self.num_comps, mc))

            cpore0_tmp = np.zeros((self.num_comps, nc, mc))#self.__altshape)
            ydot_tmp = np.zeros(self.__y0shape)
            
            def diffun(t, y0):
                nonlocal aau
                nonlocal aau2
                nonlocal cpore0_tmp
                nonlocal ydot_tmp
                idx = int(np.floor(t/tconv)) # assumes daily index provided
                extra = t/tconv - idx

                y0tmp = y0.reshape(self.__y0shape)
                ydot = ydot_tmp.copy()
                
                ### NEW Way with no loops
                ydot2 = ydot_tmp.copy() 
                
                # calculate time varying values
                cinfl = interp(cinfA[idx:idx+2], extra)
                foul_fac = interp(foulFA[idx:idx+2], extra).reshape(ThreeDSize)
                facAv = interp(facA[idx:idx+2], extra).reshape(ThreeDSize)
                
                # trying to eliminate for loops
                z2 = np.multiply(ym_vA, y0tmp[:,:nc,:mc])
                z2[z2<0] = 0.
                qte2 = np.tile(z2.sum(axis=0), ThreeDSize)  #total mass on carbon
                qte2[qte2==0.] = 1e-30 ### set mass to very small number
                yt0_2 = np.multiply(xni_vA, z2)
                yt0_2 = np.tile(yt0_2.sum(axis=0), ThreeDSize) #liquid phase in equilibrium#yt0_c
                z2 = np.divide(z2, qte2) 
                z2[qte2 <= 0.] = 0. #z_c
                
                q0_2 = np.multiply(xn_ymA, yt0_2) 
 
                cpore2 = np.power(q0_2, xni_vA) 
                cpore2 = np.multiply(z2, cpore2)
                cpore2 = np.multiply(foul_fac, cpore2)
                cpore2[qte2 <= 0.] = 0.
                cpore2[yt0_2 <= 0.] = 0. #yt0_c
                cpore2[xni_vA*np.asarray(q0_2).astype(np.float64) < -20] = 0.
                                        
                cpore_tmp2 = cpore2[:,nc-1]              
                cpore_tmp2[cpore_tmp2 < 0.] = 0.
                
                cbs2 = np.subtract(y0tmp[:,nc], cpore_tmp2) #reshape vs make outside loop
                cbs2 = np.multiply(cbs2, stdv_v.reshape(TwoDSize))
                cbs2[:,0] = 0.
                
                bb2 = np.matmul(bedp_mult[:,:-1], cpore2)
                bb2 = np.multiply(bb2, facAv)
                bb2part2 = np.matmul(beds_mult[:,:-1], y0tmp[:,:nc,:])
                bb2 = np.add(bb2, bb2part2)
                           
                ww2 = np.matmul(wr_A[:,:nd].reshape((self.num_comps,1,nd)),
                                bb2)
                
                ydot2[:,:nd,1:] = bb2[:,:,1:]
                
                num = (cinfl - cpore2[:,nc-1,0]).reshape(TwoDSize)
                num = np.multiply(num, stdv_dgI)
                num = np.subtract(num, ww2[:,0,:1])
                den = wr_A[:,nc-1].reshape(TwoDSize)
                ydot2[:, nc-1, 0] = np.divide(num,den).flatten()
                
                tmpnc1 = np.multiply(cbs2[:,1:], dgI_v)
                tmpnc1 = np.subtract(tmpnc1, 
                                     ww2[:,:,1:].reshape(self.num_comps, mc-1))
                ydot2[:, nc-1, 1:] = np.divide(tmpnc1, 
                                               wr_A[:,nc-1].reshape(TwoDSize))
                
                y0tmpV = y0tmp[:,-1,1:].reshape((self.num_comps, mc-1, 1))
                tmpaau2 = np.matmul(az_A[:,1:,1:], y0tmpV)
                aau2[:, 1:] = tmpaau2.reshape(self.num_comps, mc-1)
                
                ydotfinal = np.multiply(az_A[:,:,0],
                                        cinfl.reshape(TwoDSize))
                ydotfinal = np.add(ydotfinal, aau2)
                ydotfinal = np.multiply(np.negative(dg1), ydotfinal)
                ydotpart2 = np.multiply(3, cbs2)
                ydotfinal = np.subtract(ydotfinal, ydotpart2)
                ydot2[:,-1,1:] = ydotfinal[:,1:]
                
                ydot = ydot2.flatten()

                return ydot
 
            #Finds the locations of effluent data from results array
            effluent_locator = [(nc+1)*mc*(i+1)-1\
                                for i in range(self.num_comps)]
            
            tstart = 0 #consider moving this to self.tstart? also consider tstep
            y = solve_ivp(diffun,\
                            (tstart*tconv*t_mult, ttol),\
                            y0, \
                            method=self.solver,\
                            jac_sparsity=self.jac_sparse,\
                            # max_step=tstep/3.,\
                            )
            
            # defines interpolating function of predicted effluent
            time_index = y.t/tconv/t_mult # creates time_array and converts back to dimensions 

            cp_tmp = y.y[effluent_locator].transpose()
            cp_tmp[cp_tmp < 0.] = 0.#sets negative values to 0.
            
            cp_df = pd.DataFrame(cp_tmp, index=time_index, columns=compound_list)
            cp_df = cp_df * cb0 * mw / self.mass_mul
        
            if self.plot_output:
                cp_df.plot.line()
                
                fileID = str(np.random.randint(1e9))
                plt.savefig('fig_'+fileID+'.png', dpi=300)
                plt.close()
                
            cp = {}
            for i in range(self.num_comps):
                name = compound_list[i]
                cp[name] = interp1d(y.t/tconv/t_mult, 
                                    cp_df[name], 
                                    fill_value='extrapolate')
            return cp
        
        best_fit = run()
        return best_fit
    
##### end run_psdm()    

    def run_psdm_kfit(self, compound):
        ## now uses run_psdm()
        idx=pd.IndexSlice     
        
        mp.freeze_support()
        k = 1
        xn = 1 
        ssqs = 1
        results = 1
        
        # print(self.data_df)

        self.num_comps = 1
        self.jac_sparse = spar_Jac(self.num_comps, self.nc, self.nz, self.ne)
        self.__y0shape = (self.num_comps, self.nc+1, self.mc)
        self.__altshape = (self.num_comps, self.nc, self.mc)
        
        self.compounds = [compound]
            
        ## get influent/effluent data for single species
        self.data_df = self.data_bup.transpose().loc[idx[:, compound], :].transpose()
        
        # self.data_df = filter_compounds(self.data_df, [compound], self.carbon, self.influent)
        f_eff = interp1d(self.data_df.index, self.data_df[self.carbon, compound], fill_value='extrapolate')
        
        ssq_xs = np.arange(self.k_data[compound]['brk']) ## only consider through breakthrough, old: np.max(self.data_df.index))
        
        ssqs = pd.DataFrame(index=self.test_range, columns=self.xn_range)
        
        if len(self.test_range) == 1 and len(self.xn_range) == 1:
            self.optimize_flag = False
        
        ## testing optimized loop 
        if self.optimize_flag:
            ### should this if be removed, run_all and run_all_smart basically handle this now
            
            if compound not in self.k_by_xn_factor.keys():
                ### hopefully not really used
                k_mult = {}
                for xns in self.xn_range:
                    k_mult[xns] = recalc_k(self.k_data[compound]['K'], 
                                           self.k_data[compound]['q'],
                                           self.k_data[compound]['1/n'],
                                           xns)
            else:
                ### use precalculated k_multipliers
                k_mult = {i: self.k_by_xn_factor[compound](i) for i in self.xn_range}
            
            for i in self.test_range:
                for xn in self.xn_range:
                    self.k_data.loc['1/n', compound] = xn
                    self.k_data.loc['K', compound] = k_mult[xn] * i
                    
                    results = self.run_psdm()
                    
                    ssq = ((results[compound](ssq_xs) - f_eff(ssq_xs))**2).sum() ## returns sum of squares for assumed daily data with linear interpolated effluent
                    ssqs.loc[i, xn] =ssq
           
        else:
            # print(self.k_data_input_type)
            if self.k_data_input_type == list:
                if len(self.test_range) == 1 and len(self.xn_range) == 1:
                    ### This means nothing was input and input occured through test_range and xn_range
                    self.k_data.loc['1/n', compound] = self.xn_range[0]
                    self.k_data.loc['K', compound] = self.test_range[0]
                    
                    results = self.run_psdm()
                    
                    ssq = ((results[compound](ssq_xs) - f_eff(ssq_xs))**2).sum() ## returns sum of squares for assumed daily data with linear interpolated effluent
                    ssqs.loc[self.test_range[0], self.xn_range[0]] = ssq
                    
                    results = pd.DataFrame(results[compound].y, index=results[compound].x,
                                           columns=[compound])
                
                else:
                    ### try without changing K & 1/n?
                    results = self.run_psdm()
                   
                    ssq = ((results[compound](ssq_xs) - f_eff(ssq_xs))**2).sum() ## returns sum of squares for assumed daily data with linear interpolated effluent
                   
                    k = self.k_data[compound]['K']
                    xn = self.k_data[compound]['1/n']
                    ssqs = pd.DataFrame(ssq, index=[xn],
                                            columns=[k])
                   
                    results = pd.DataFrame(results[compound].y, index=results[compound].x,
                                           columns=[compound])
                   
            else:
                ## This should mean that a k_data DataFrame was input, can just run
                results = self.run_psdm()
                
                ssq = ((results[compound](ssq_xs) - f_eff(ssq_xs))**2).sum() ## returns sum of squares for assumed daily data with linear interpolated effluent
                
                k = self.k_data[compound]['K']
                xn = self.k_data[compound]['1/n']
                ssqs = pd.DataFrame(ssq, index=[xn],
                                    columns=[k])
                
                results = pd.DataFrame(results[compound].y, index=results[compound].x,
                                       columns=[compound])
                
        
        ### END, cleanup
        self.__reset_column_values()

        return compound, k, xn, ssqs, results
        ## end run_psdm_kfit()
        
    def run_psdm_dsfit(self, compound):
        idx=pd.IndexSlice     
        ## replacement of run_psdm_dsfit() that uses run_psdm() ### TODO: Need to test
        mp.freeze_support()
 
        ssqs = 1
        results = 1

        ## always makes it a single case
        self.num_comps = 1
        self.jac_sparse = spar_Jac(self.num_comps, self.nc, self.nz, self.ne)
        self.__y0shape = (self.num_comps, self.nc+1, self.mc)
        self.__altshape = (self.num_comps, self.nc, self.mc)
        
        self.compounds = [compound]
            
        ## get influent/effluent data for single species
        self.data_df = self.data_bup.transpose().loc[idx[:, compound], :].transpose()
        
        # self.data_df = filter_compounds(self.data_df, [compound], self.carbon, self.influent)
        f_eff = interp1d(self.data_df.index, self.data_df[self.carbon, compound], fill_value='extrapolate')

        mol_vol = self.k_data[compound]['MolarVol']
        mw = self.k_data[compound]['MW']
        molar_k = self.k_data[compound]['K'] / mw / ((1. / mw) ** self.k_data[compound]['1/n'])

        cb0 = self.data_df[self.influent, compound][0] * self.mass_mul / mw

        difl = 13.26e-5/(((self.vw * 100.)**1.14)*(mol_vol**0.589)) #vb
        ds_base = self.epor * difl * cb0 * self.psdfr / (1e3 * self.rhop * cb0**self.k_data[compound]['1/n'])

        # tests = self.test_range
        tests = np.linspace(1e-10, 1, num=30)

        ssqs = pd.Series(index=tests) ## only consider test range in ds direction
        best_ssq = 1e20 ## big initial value
        best_factor = 1

        ### Ds fit only used for fitting, so ignore self.optimize_flag
        for factor in tests: ### need to swap out
            self.mass_transfer.loc['ds', compound] = factor * ds_base
            ### generally test_range should be something like np.linspace(1e-4, 1, num=20), where Ds is likely to be lower than predicted
            results = self.run_psdm()
                    
            ssq = ((results[compound](ssq_xs) - f_eff(ssq_xs))**2).sum() ## returns sum of squares for assumed daily data with linear interpolated effluent
            ssqs.loc[factor] =ssq

            if ssq < best_ssq:
                best_ssq = ssq * 1
                best_factor = factor * 1
                results_out = pd.DataFrame(results[compound].y, index=results[compound].x,
                                       columns=[compound])



        # print(self.mass_transfer)
                
        
        ### END, cleanup
        self.__reset_column_values()

        return compound, best_val_ds, ssqs, results_out, ds_base
        ## end run_psdm_dsfit

    ### Begin model_uncertainty()
    def model_uncertainty(self, single=True, capacity=10, k='None', qn='None', c0='None', mass='None', flrt='None', ds='None', dp='None', kf='None'):
        idx=pd.IndexSlice 
        
        self.__set_backups()

        ## create results dataframe: self.results stores baseline results, without considered uncertainty
        self.results = pd.DataFrame(index=np.linspace(0, self.duration, num=2000), columns=self.compounds_bup)
        ### create uncertainty_results dataframe
        multi_idx = [(i, j) for i in ['upper', 'lower'] for j in self.compounds_bup]
        midx = pd.MultiIndex.from_tuples(multi_idx)
        self.uncertainty_results = pd.DataFrame(index=self.results.index, columns=midx)
        
        calced_mass_transfer = self.mass_transfer.copy()
        
        #### run PSDM model from inputs - no changes
        if single:    ## run compounds as single compounds, default  
            self.num_comps = 1
            self.jac_sparse = spar_Jac(self.num_comps, self.nc, self.nz, self.ne)
            self.__y0shape = (self.num_comps, self.nc+1, self.mc)
            self.__altshape = (self.num_comps, self.nc, self.mc)
            
            for comp in self.compounds_bup: ## run over all compounds
                self.compounds = [comp]
                
                ## get influent/effluent data for single species
                self.data_df = self.data_bup.transpose().loc[idx[:, comp], :].transpose()
                
                ## set k_data
                self.k_data = pd.DataFrame(self.k_data_bup[comp].values, index=self.k_data_bup.index, columns=self.compounds) 

                temp_results = self.run_psdm() ## run simulation
                self.results[comp] = temp_results[comp](self.results.index) ## store results in dataframe
                
                calced_mass_transfer[comp] = self.mass_transfer_data[comp]
                #print(self.mass_transfer_data[comp])
        else: ## run as multi-competitive
        ## TODO: Still need to test this, but it should work...
            temp_results = self.run_psdm() ## run simulation, returns dictionary of interpolating functions
                
            for comp in temp_results.keys(): ## store results in dataframe
                self.results[comp] = temp_results[comp](self.results.index)
        
            calced_mass_transfer = self.mass_transfer_data.copy()
        
        
        ## begin uncertainty loop #######################################################################
        
        ### Create copy of base results to start uncertainty_results
        for comp in self.compounds_bup:
            # print(comp)
            self.uncertainty_results.loc[self.results.index, ('upper', comp)] = self.results[comp] * 1
            self.uncertainty_results.loc[self.results.index, ('lower', comp)] = self.results[comp] * 1
       
        test_uncertainty = [] ## initialize list, will be list of dictionaries 
        
        ### Inputs assume % reported as 10 for 10%, not 0.1
        compare_bounds = False
        if capacity != 'None':
            test_uncertainty.append({'k': 1 + capacity/100, '1/n': 1 - capacity/100}) ## Freundlich K and 1/n work in opposite directions
            test_uncertainty.append({'k': 1 - capacity/100, '1/n': 1 + capacity/100})
            compare_bounds = True
        
        other_inputs = {'k': k, 'qn': qn, 'mass': mass, 'ds': ds, 'dp': dp, 'kf': kf, 'flrt': flrt}
        for key, value in other_inputs.items():
            
            if value != 'None':
                test_uncertainty.append({key: 1 + value/100})
                test_uncertainty.append({key: 1 - value/100})
                compare_bounds = True
        

        for test in test_uncertainty:
            self.__reset_column_values()
            k_data_loop = self.k_data.copy()
            self.mass_transfer = calced_mass_transfer.copy() ## this always fixes mass transfer to original system flowrates... TODO: is that a problem?
            data = self.data_df
            
            for key in test.keys(): ## this should be a dictionary... 
                
                if key == 'k':
                    k_data_loop.loc['K'] *= test[key]
                
                if key == '1/n':
                    k_data_loop.loc['1/n'] *= test[key]
                
                if key == 'flrt': ## adjust flowrate
                    self.flrt *= test[key]
                    
                if key == 'ds':  ## adjust surface diffusion coefficient
                    self.mass_transfer.loc['ds'] *= test[key]
                    ## This is only really relevant if fouling is not used. Otherwise, Ds becomes negligibly small
                
                if key == 'dp':  ## adjust pore diffusion coefficient
                    self.mass_transfer.loc['dp'] *= test[key]
                
                if key == 'kf': ## adjust film transfer coefficient
                    self.mass_transfer.loc['kf'] *= test[key]
            
                if key == 'mass': ## adjust mass of media
                    self.wt *= test[key]
                    
                    volume = np.pi/4 * (self.diam**2) * self.L
                    
                    if self.wt/volume > self.rhof:
                        self.L = self.wt / (np.pi/4 * self.diam**2 * self.rhof) ## may cause shift to EBCT, scaling will always expect base case as EBCT adjustment
                        ## resize column to make sure apparent density doesn't exceed rhof
                        
                if key == 'C': ## adjust influent concentration
                    data[self.influent] *= test[key] 
                    self.data_df = data.copy()  # overwritten if "single=True", but used by competitive if not.
                    self.k_data.loc['AveC'] *= test[key] ## adjust average C for molar K adjustment
                
                if key == 'qn':  ### adjust capicity (q's) impact on 1/n
                    
                    qs = self.k_data.loc['q'] * test[key]
                    cs = self.k_data.loc['AveC']
                    ks = self.k_data.loc['K']
                    
                    self.k_data.loc['1/n'] = np.log(qs / ks) / np.log(cs * self.mass_mul)
        
            ### calculate updated values if uncertainty has changed these values
            #calculate other fixed values
            self.area = np.pi*(self.diam**2)/4.
            self.bedvol = self.area * self.L
            self.ebed = 1. - self.wt/(self.bedvol*self.rhop)
            self.tau = self.bedvol * self.ebed * 60./self.flrt
            self.sf = 0.245423867471 * self.flrt/self.area # gpm/ft**2
            self.vs = self.flrt/(60.*self.area)
            self.re = (2.*self.rad*self.vs*self.dw)/(self.ebed*self.vw)
            
            #calculate Empty Bed Contact Time (EBCT)
            self.ebct = self.area*self.L/self.flrt 
            
            ## run uncertainty
            if single:    ## run compounds as single compounds, default  
                
                self.num_comps = 1
                self.jac_sparse = spar_Jac(self.num_comps, self.nc, self.nz, self.ne)
                self.__y0shape = (self.num_comps, self.nc+1, self.mc)
                self.__altshape = (self.num_comps, self.nc, self.mc)
                
                for comp in self.compounds_bup: ## run over all compounds
                    self.compounds = [comp]
                    
                    ## get influent/effluent data for single species
                    self.data_df = data.transpose().loc[idx[:, comp], :].transpose()
                    
                    ## set k_data
                    self.k_data = pd.DataFrame(k_data_loop[comp].values, index=k_data_loop.index, columns=self.compounds) 
                    
                    temp_results = self.run_psdm() ## run simulation
                    
                    ### figure out if this should change the uncertainty bounds
                    if compare_bounds:
                        ## only check if uncertainty is modeled
                        lower = np.minimum(self.uncertainty_results['lower'][comp].values, temp_results[comp](self.results.index))
                        upper = np.maximum(self.uncertainty_results['upper'][comp].values, temp_results[comp](self.results.index))
                        
                        self.uncertainty_results.loc[idx[:, ('lower', comp)]] = lower * 1
                        self.uncertainty_results.loc[idx[:, ('upper', comp)]] = upper * 1 ## store results in dataframe
                    
            else: ## run as multi-competitive
            ## TODO: Still need to test this, but it should work...
            
                self.k_data = k_data_loop.copy()
                temp_results = self.run_psdm() ## run simulation, returns dictionary of interpolating functions
                    
                for comp in temp_results.keys(): ## store results in dataframe
                    self.results[comp] = temp_results[comp](self.results.index)
                    
                    if compare_bounds:
                        ## only check if uncertainty is modeled
                        lower = np.minimum(self.uncertainty_results['lower'][comp].values, temp_results[comp](self.results.index))
                        upper = np.maximum(self.uncertainty_results['upper'][comp].values, temp_results[comp](self.results.index))
                        
                        self.uncertainty_results['lower'][comp] = lower
                        self.uncertainty_results['upper'][comp] = upper ## store results in dataframe
            
        ### end uncertainty loop ##########################################################################

        
        ## Final reset to ensure things are back to original state
        self.__reset_column_values()

        ### no value returned, saved as self.results and self.uncertainty_results
        

### end model_uncertainty()





    
# =============================================================================
#     NEW BELOW, may delete
# =============================================================================
      
    
    def run_all_MP(self, plot=False, save_file=False):
        '''
        Parameters
        ----------
        plot : BOOL, optional
            Are plots generated during this function. The default is False.
        save_file : BOOL, optional
            Are results files generated during this function. The default is True.
        optimize : string, optional
            'brute' or 'staged' are acceptable options. The default is 'staged'.
        init_grid : INT, optional
            Number of grid points for staged optimizer. The default is 5.
        init_loop : INT, optional
            Number of refinement loops for staged optimizer. The default is 3.

        Current Behavior
        ----------------
        All files will be overwritten. No file checking is currently performed.
        
        Returns
        -------
        None. Optimized k & 1/n values are saved to Class object, 
        and can be saved in script that calls this function.

        '''
        
        
        #forces single run to handle optimizing in this funciton, not run_psdm_kfit
        opt_flg = self.optimize_flag
        orig_test_range = self.test_range * 1.
        orig_xn_range = self.xn_range * 1.
        self.optimize_flag = False
        
        num_processes = self.processes - 2
        
        k_range = self.test_range
        xn_range = self.xn_range


        for compound in self.compounds:
            print(compound, ' running')

            k_val = self.k_data[compound]['K']
            q_val = self.k_data[compound]['q']
            
            inf = self.data_df[self.influent][compound]
            eff = self.data_df[self.carbon][compound]
            
            k_mult = {}
            models = []
            
            
            for xns in xn_range:
                k_mult[xns] = recalc_k(k_val, q_val, self.xn, xns)
            
            for k in k_range:
                for xn in xn_range:
                    models.append([deepcopy(self), k * k_mult[xn], xn, 
                                   compound, k])

            # processes = []
            # for row in models:
            #     p = mp.Process(target=run_MP_helper, args=(*row, ))
            #     processes.append(p)
            
            # result_list = [x.start() for x in processes]
            # print('Here')
            # for item in result_list:
            #     item.join()
            
            pool = mp.Pool(processes=num_processes)
            retval = pool.starmap(run_MP_helper, models)
            pool.close()
            pool.join()
                    
            data_pd = pd.DataFrame(retval, columns=['K2','1/n','ssq','compound','K'])
            ssqs = pd.pivot_table(data_pd,values=['ssq'], 
                                             index=['K'],
                                             columns=['1/n'],
                                             aggfunc=np.max,
                                             )['ssq']
            
            if ~np.isclose(ssqs.min().min(), ssqs.max().max()):# np.floor(ssqs.values[0][0]):
                min_val = find_minimum_df(ssqs)
                best_val_xn = min_val.columns[0]
                best_val_k = min_val.index[0] * k_mult[best_val_xn]
                min_val = min_val.values[0][0]
                
                self.xn_range = np.array([best_val_xn])
                self.test_range = np.array([best_val_k])
                _, _, _, _, md = self.run_psdm_kfit(compound)
                
                if plot:
                    plt.figure()
                    plt.plot(inf.index, inf.values, marker='x', label='influent')
                    plt.plot(eff.index, eff.values, marker='.',
                              markerfacecolor='None',label='effluent')
                    #plot model results
                    plt.plot(md.index, 
                              md.values,
                              label=f"K: {best_val_k:.2f}, 1/n: {best_val_xn:.3f}")
    
                    
                    plt.legend()
                    plt.title(f"{compound} - {self.carbon}")
                    plt.savefig(f"{self.carbon}_{compound}.png",dpi=300)
                    plt.close()
                    
                    plt.figure()
                    ssqs[ssqs>np.percentile(ssqs, 25)] = np.percentile(ssqs, 25)
                    plt.contourf(ssqs.columns, ssqs.index, ssqs.values)
                    plt.title(f"{compound} - {self.carbon}")
                    plt.savefig(f"ssq_{self.carbon}_{compound}.png",dpi=300)
                    plt.close()
                
                if save_file:
                    with pd.ExcelWriter(f"ssq_{self.carbon}-{compound}.xlsx") as writer:
                        ssqs.to_excel(writer, 'Sheet1')
                    
                    with pd.ExcelWriter(self.project_name+'_'+self.carbon + '-' + \
                                            compound+'.xlsx') as writer:
                
                        md.to_excel(writer, 'model_fit')
                        
                        inf.to_excel(writer, 'influent')
                        eff.to_excel(writer, 'effluent')
                        
                        data_tmp = pd.Series([self.re, best_val_k, best_val_xn,\
                                      min_val, self.ebct, self.sf], \
                                      index=['Re','K','1/n','ssq','ebct','sf'])
                        data_tmp.to_excel(writer, 'parameters')
                
                self.k_data[compound]['K'] = best_val_k * 1
                self.k_data[compound]['1/n'] = best_val_xn * 1

        self.xn_range = orig_xn_range
        self.test_range = orig_test_range
        self.optimize_flag = opt_flg #resets to original value

    
    
# =============================================================================
# END OF PSDM CLASS
# =============================================================================
