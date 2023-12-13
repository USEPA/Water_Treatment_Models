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
import scipy as sp
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import interp1d
from scipy.stats import linregress
import multiprocessing as mp
import time as ti #time as a variable is used in code, so ti is used

#Read in all associated PSDM functions
from PSDM_functions import min_per_day, lpg, spar_Jac, foul_params, kf_calc
from PSDM_functions import find_minimum_df, tortuosity, calc_solver_matrix
from PSDM_functions import density, viscosity, recalc_k, generate_grid
from PSDM_functions import interp 
from PSDM_functions import process_input_data, process_input_file

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
        self.temp = kw.get('temp',20)
        self.vw = viscosity(self.temp)
        self.dw = density(self.temp)
        
        #define unit conversions (time)
        if 'time' not in column_data.index:
            self.time_type = kw.get('time_type', 'days') 
            if self.time_type == 'days': #base units in minutes
                self.t_mult = min_per_day
            elif self.time_type == 'hours':
                self.t_mult = 60.
            else:
                self.t_mult = 1.
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
                self.flow_mult = 1.       
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
            else:
                print('conc_type is invalid, ug/ng are valid options')
                self.mass_mul = 1.        
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

        #calculate initial values 
        #might need to add more here
        k_data = kw.get('k_data', [])
        if len(k_data) == 0:
            self.k_data = pd.DataFrame(index=['K','1/n', 'q', 'brk','AveC'], \
                                        columns=self.compounds)
            for comp in self.compounds:
                k, q, classifier, brk, aveC = self.__calculate_capacity(comp)
                self.k_data[comp]=np.array([k, self.xn, q, brk, aveC])
        else:
            self.k_data = k_data
            if self.brk_type=='force':
                tmp_assume = self.brk_df[(self.brk_df['carbon']==self.carbon)&\
                                         (self.brk_df['breakday']!=3000)]
                impacted = tmp_assume['compound'].values
                for comp in self.compounds:
                    if comp in impacted:
                        k, q, classifier, brk, aveC = self.__calculate_capacity(comp)
                        self.k_data[comp]=np.array([k, self.xn, q, brk, aveC])
            
        self.k_data_orig = self.k_data
        
        self.water_type = kw.get('water_type','Organic Free')
        self.chem_type = kw.get('chem_type', 'halogenated alkenes')
        
        self.test_range = kw.get('test_range', np.linspace(1, 5, 41))
        self.xn_range = kw.get('xn_range', np.arange(0.20, 0.95, 0.05))
        
        #handling for multiprocessing
        self.processes = kw.get('mp', mp.cpu_count())
        
        self.optimize_flag = kw.get('optimize',True)
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
        self.fouling_dict= self.__fouled_k_new(self.k_data.loc['K'], self.time_vals)
        
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
    
    def __fouled_k_new(self, k_data, t):
        #works on unconverted time
        #only works on multiplier
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
        flow = self.flrt * self.flow_mult
        breakthrough_code = 'none'
        carbon_mass = self.wt
        infl = self.data_df[self.influent][compound]
        infl[infl == 0] = 1e-3  ### prevents divide by zero error
        effl = self.data_df[self.carbon][compound]
        
        breakthrough_time = self.duration
        k = 0.
        q_meas = 0.
        xdata = self.xdata
        aveC = infl.mean()              #calculates average influent concentration
        
        if self.brk_type == 'calc':
            donothing = False #set trigger for later process to False
            if infl.sum() == 0:
                print('No mass in influent for '+compound)
                donothing = True
            if effl.sum() == 0:
                print('Insufficient data for breakthrough for '+compound)
                donothing = False    
            
            diff = (infl-effl)#/infl
            perc_diff = (infl-effl)/infl
        
            done = False
            if not donothing:
                yte = effl.values.astype('float64')
                yti = infl.values.astype('float64')
                
                #check if there are already zeros (breakthrough point exact)
                tmp_time = perc_diff[perc_diff==0.].index
                if len(tmp_time) == 1:
                    breakthrough_time = tmp_time.values[0]
                    breakthrough_code = 'breakthrough'
                    done = True
                elif len(tmp_time) > 1:
                    breakthrough_time = tmp_time.min()
                    breakthrough_code = 'breakthrough'
                    done = True
                #check if there are transitions to negative (breakthrough in data, calculatable)
                if not done:
                    tmp_time = diff[diff<0.].index.values
                    if len(tmp_time) > 0:
                        upper = diff[diff >= 0]#.iloc[-1]
                        upper_x = upper.index[-1]
                        upper_y = upper.iloc[-1]
                        lower = diff[diff < 0]#.iloc[0]
                        lower_x = lower.index[0]
                        lower_y = lower.iloc[-1]
                        
                        # print(upper_x, upper_y, lower_x, lower_y)
                        slope = (lower_y - upper_y) / (lower_x - upper_x)
                        intercept = upper_y - slope * lower_x
                        breakthrough_time = -intercept/slope
                        
                        breakthrough_code = 'breakthrough'
                        done = True
                
                #check to see if there is possible breakthrough within 15% of influent
                #and minimal change in derivative
                if not done:
                    tmp_time = perc_diff[perc_diff < 0.1].index.values #
                    if len(tmp_time) > 0:
                        test = np.min(tmp_time)
                        if test != np.max(xdata):
                            tmp_diff = perc_diff[perc_diff.index>=test]
                            if np.max(tmp_diff.values) < 0.15 and np.min(tmp_diff-tmp_diff.values[0])>=0.:
                                breakthrough_code = 'implied'
                                breakthrough_time = test
                                done = True
                
                ### check for correlational agreement with values, to imply breakthrough
                if not done:
                    length = len(yte)
                    corrv = np.array([np.corrcoef(yte[i:i+3],yti[i:i+3])[0,1] \
                                    for i in range(length-2)])
                    corrv[np.isnan(corrv)] = 0.        
                    if corrv[-1] >= 0.95 and yte[-1] > 0.:
                        for i in range(len(corrv)-1,0,-1):
                            if corrv[i] >= 0.95:
                                breakthrough_time = xdata[i]
                                breakthrough_code = 'implied'
                                done = True
                            elif corrv[i] < 0.95:
                                break
        

                #try to estimate the breakthrough
                if not done:                
                    nZc = np.count_nonzero(yte)     #used to determine linear regression

                    if nZc < len(yte):
                        if yte[-(nZc+1)] != 0:
                            nZc += 1 ## adds one more to non-zero value, if original nZc does not equal 0.
                        x = xdata[-(nZc+1):].astype('float64')
                        y = yte[-(nZc+1):].astype('float64')
                    else:
                        x = xdata[-(nZc):].astype('float64')
                        y = yte[-(nZc):].astype('float64')
                    

                    
                    m, b, *extra = linregress(x, y)
                    
                    if m > 0:
                        intersection = (infl.mean() - b)/m
                    else:
                        intersection = self.duration
                    
                    breakthrough_time = intersection
                    breakthrough_code = 'estimated'
                    
                    yti = np.append(yti, aveC)
                    yte = np.append(yte, aveC)
                    xdata = np.append(xdata, breakthrough_time)
           
            
                f_infl = interp1d(xdata, yti, fill_value='extrapolate')
               
                if breakthrough_time <= np.max(xdata):
                    xdata_trunc = xdata[xdata <= breakthrough_time]
                    num_vals = len(xdata_trunc)
                    f_effl = interp1d(xdata[:num_vals+1], yte[:num_vals+1], fill_value='extrapolate')
                else:
                    f_effl = interp1d(xdata, yte, fill_value='extrapolate')
                
                try:
                    int_infl = quad(f_infl, 0, breakthrough_time, points=xdata)[0]
                    int_effl = quad(f_effl, 0, breakthrough_time, points=xdata)[0]
                except Exception as e:
                    print(e)
                    int_infl = quad(f_infl, 0, breakthrough_time)[0]
                    int_effl = quad(f_effl, 0, breakthrough_time)[0]
                qtmp = flow * self.t_mult * (int_infl - int_effl) * self.mass_mul
                q_meas = qtmp/carbon_mass #ug/g
                aveC = int_infl/breakthrough_time # recalculate only what is inside breakthrough
                k = q_meas / ((aveC*self.mass_mul) ** self.xn) 
                
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
                qtmp = flow * self.t_mult * (int_infl - int_effl) * self.mass_mul
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
                        
                q_meas = flow * self.t_mult * qtmp * self.mass_mul /\
                         (carbon_mass) # ug/g
                k = q_meas / ((aveC * self.mass_mul) ** self.xn) 
                
            breakthrough_code = 'supplied'
            
        # returns capacity in (ug/g)(L/ug)**(1/n)
        return k, q_meas, breakthrough_code, breakthrough_time, aveC 
    
# =============================================================================
# End INIT and helper functions
# =============================================================================
    
    def run_psdm_kfit(self, compound):
        '''
        time: must be specified in minutes
        best_vals: [Dp, Ds, kf]
        flow rate assumed to be in 'ml/min' in column properties
        '''
        mp.freeze_support()
        
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
        molar_k_t = self.fouling_dict[compound] #multiplier for fouling
        
        water_type = self.water_type
        
        k_v = self.k_data[compound]['K']
        q_v = self.k_data[compound]['q']
        
        mw = self.comp_df[compound]['MW']           #molecular weight
        mol_vol = self.comp_df[compound]['MolarVol'] # Molar volume
        
        if self.time_type == 'days':
            dstep = 0.25 * min_per_day
        else:
            dstep = 15.        
        
        #set up bindings for nonlocal varaibles
        cinf = 1.
        cout_f = 1.
        tconv = 1.
        time_dim = 1.
        time_dim2 = 1.
        ttol = 1.
        tstep = 1.
        ds_v = 1.
                
        inf = self.data_df[self.influent][compound]
        eff = self.data_df[self.carbon][compound] 
        #convert cbo to molar values
        cbo = inf * self.mass_mul / mw 
        time = (inf.index * t_mult).values
        if inf.index[-1] < 10 and self.time_type == 'days':
            dstep = 15.
        elif self.time_type == 'min':
            dstep = 1.
            
        if cbo.iloc[0] == 0.:
            cb0 = 1.
        else:
            cb0 = cbo[0]
        
        cin = cbo/cb0
        
        try:
            brk = self.k_data[compound]['brk']
        except:
            brk = np.max(inf.index.values)
               
        tortu = self.tortu                             # tortuosity
        psdfr = self.psdfr                             # pore to surface diffusion ratio
        nd = nc - 1
        
        difl = 13.26e-5/(((vw * 100.)**1.14)*(mol_vol**0.589)) #vb
        sc = vw / (dw * difl)       #schmidt number
        
        #set film and pore diffusion
        multi_p = difl/(2*rad) # multiplier used for kf calculation
        if self.mass_transfer[compound]['kf'] == 0:
            kf_v = kf_calc(multi_p, self.re, sc, ebed, corr='Chern and Chien')
        else:
            kf_v = self.mass_transfer[compound]['kf']
            

        if compound == 'Test':
            kf_v = self.k_data['Test']['kf'] #will break
            
        cout = eff * self.mass_mul / mw / cb0 
        
        if self.mass_transfer[compound]['dp'] == 0.:
            dp_v = (difl/(tortu))       #porosity not used in AdDesignS appendix, removed to match
        else:
            dp_v = self.mass_transfer[compound]['dp']
        
        def run(k_val, xn):
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
            molar_k = k_val / mw / ((1. / mw) ** xn)  
            xni = 1./xn
            
            
            if self.mass_transfer[compound]['ds'] == 0.:
                ds_v = epor*difl*cb0*psdfr/(1e3*rhop*molar_k*cb0**xn)
            else:
                ds_v = self.mass_transfer[compound]['ds']
            
            
            if water_type != 'Organic Free':
                ds_v /= 1e10 #1e6

            d = ds_v/dp_v
            
            qe = molar_k * cb0**xn 
            qte = 1. * qe
            
            dgs = (rhop * qe * (1.-ebed) * 1000.)/(ebed * cb0)
            dgp = epor * (1. - ebed)/(ebed) 
            dg = dgs + dgp
            dgt = dg
            dg1 = 1. + dgt
            dgI = 1.0/dg
            edd = dgt/dg #dgt changed from dg1
            
            ym = qe/qte
            
            eds = ds_v*dgs*tau/(rad**2)
            if eds < 1e-130:
                eds = 1e-130
            edp = dp_v*dgp*tau/(rad**2)
            # from orthog(n)
            beds = (eds + d*edp) * edd * br[:-1]
            bedp = edp * (1. - d) * edd * br[:-1]
            #depends on kf
            st = kf_v * (1. -ebed) * tau/(ebed*rad)
            stdv = st * dgt     # dgt changed from dg1
            
            #convert to dimensionless
            tconv = 60./(tau*dg1)
            self.tconv = tconv
            tstep = dstep * tconv
            ttol = time[-1] * tconv
            time_dim = time * tconv
            
            numb = int(brk*3 + 1)
            time_dim2 = np.linspace(0., brk * tconv * t_mult,\
                                    num=numb, endpoint=True) #increase the number of sites to check for ssq
            #set up time based influent data
            cinf = interp1d(time, cin.values, fill_value='extrapolate') 
            # cinf = interp1d(time_dim, cin.values, fill_value = 'extrapolate') 
            cout_f = interp1d(time_dim, cout.values, fill_value='extrapolate')
            
            #initialize storage arrays/matrices
            n = (nc+1)*mc
            y0 = np.zeros(n)
            
            #new array methods, create all arrays needed by diffun
            time_temp = np.arange(0, time[-1] + self.t_mult * 10, 1)
            cinfA = cinf(time_temp)
            if water_type != 'Organic Free':
                tortu = tortuosity(time_temp)
            else:
                tortu = np.ones(len(time_temp)) * self.tortu ## maybe?
            facA = (1./tortu - d)/(1. - d)
            foulFA = (1./molar_k_t(time_temp))**xni
            
            ydot_tmp = np.zeros((nc+1, mc))
            self.ydot = ydot_tmp * 1.
            
            def diffun(t, y0):
                nonlocal aau
                nonlocal ydot_tmp

                y0tmp = y0.reshape(ydot_tmp.shape)
                ydot = ydot_tmp.copy()
                
                idx = int(np.floor(t/tconv)) # assumes daily index provided
                extra = t/tconv - idx
                
                # #defines the influent concentration at time t
                cinfl = interp(cinfA[idx: idx+2], extra) # 
                
                z = ym * y0tmp[:nc, :mc] #updated ym should always be 1 for single comp.
                qte = z
                yt0 = xni * z
                
                z_c = z/qte
                z_c[qte<=0.] = 0.
                # z[qte>0.] = z_c[qte>0.] # should be 1 for single component.
                q0 = yt0 * xn/ym

                cpore = z_c * q0**xni * interp(foulFA[idx:idx+2], extra) 
                cpore[np.logical_or.reduce((qte<=0.,\
                                            yt0<=0,\
                                            xni*np.log10(q0)<-20,\
                                            ))] = 0.
                cpore_tmp = cpore[nc-1]
                cpore_tmp[cpore_tmp < 0.] = 0.
                cbs = stdv*(y0tmp[nc]-cpore_tmp)
                cbs[0] = 0. 
                
                bb = interp(facA[idx:idx+2], extra)*np.dot(bedp, cpore) +\
                      np.dot(beds, y0tmp[:nc, :])
                
                ww = wr[:nd]@bb

                ydot[:nd,1:] = bb[:, 1:]
            
                ydot[nc-1][0] = (stdv*dgI*(cinfl - cpore[nc-1][0]) - ww[0])/\
                                wr[nc-1] #iii
                ydot[nc-1][1:] = (cbs[1:]*dgI - ww[1:])/wr[nc-1]
                
                aau[1:] = (np.dot(az[1:,1:], y0tmp[-1, 1:]))
                
                ydot[-1,1:] = (-dgt*(az[:,0]*cinfl + aau) - 3.* cbs)[1:]  #dgt was changed from dg1  
                ydot = ydot.reshape((nc+1)*(mc))
                return ydot
            
            try:
                y = solve_ivp(diffun,\
                                (0, ttol),\
                                y0, \
                                method=self.solver,\
                                jac_sparsity=self.jac_sparse_single,\
                                max_step=tstep/3,\
                                )
                # defines interpolating function of predicted effluent
                cp_tmp = y.y[-1]
                cp_tmp[cp_tmp < 0.] = 0.#sets negative values to 0.
                cp = interp1d(y.t, cp_tmp, fill_value='extrapolate') 
                self.ydot = y.y * cb0 * mw / self.mass_mul
                self.yt = y.t / tconv / t_mult
            except Exception as e:
                print(f"Error produced: {compound} - {e}")
                t_temp = np.linspace(0, ttol, 20)
                cp_tmp = np.ones(20) # need a better error position
                cp = interp1d(t_temp, cp_tmp, fill_value='extrapolate')
                ### below 2 lines are still causing issues with lead_lag code
                self.yt = t_temp / tconv / t_mult
                self.ydot = np.zeros(((nc+1)*(mc), len(self.yt))) ### should at least make them the same shape as expected
            
            ssq = ((cout_f(time_dim2)-cp(time_dim2))**2).sum()
            return cp, ssq
        
        def run_fit(k_val, xn):
            cp, ssq = run(k_val, xn)
            return ssq
        
        if self.optimize_flag:
            k_mult = {}
            for xns in self.xn_range:
                k_mult[xns] = recalc_k(k_v, q_v, self.xn, xns)
                
            ssqs = pd.DataFrame([[run_fit(i*k_mult[j],j)\
                                  for j in self.xn_range] \
                                  for i in self.test_range], \
                                  index=self.test_range, \
                                  columns=self.xn_range)
            min_val = find_minimum_df(ssqs)
            best_val_xn = min_val.columns[0]
            best_val_k = min_val.index[0] * k_mult[best_val_xn]
            best_fit, _ = run(best_val_k, best_val_xn)
            min_val = min_val.values[0][0]
            
            with pd.ExcelWriter('ssq_'+self.carbon+'-'+compound+'.xlsx') as writer:
                ssqs.to_excel(writer, 'Sheet1')
        
        else: #assume test_range and xn_range are single values
            best_val_xn = self.xn_range[0]
            best_val_k = self.test_range[0]
            best_fit, min_val = run(best_val_k, best_val_xn)
            ssqs = pd.DataFrame(min_val, columns=[best_val_xn],\
                                index=[best_val_k])
        
        itp = np.arange(0., time[-1]/t_mult, dstep/t_mult)
        output_fit = interp1d(itp, \
                              best_fit(itp*tconv*t_mult) * cb0 *\
                              mw / self.mass_mul, \
                              fill_value='extrapolate')
        model_data = pd.DataFrame(output_fit(itp), \
                                  columns = ['data'], \
                                  index = itp)
        
        data_tmp = pd.Series([sc, self.re, difl, kf_v, best_val_k, best_val_xn,\
                     dp_v, ds_v, min_val, self.ebct, self.sf], \
                     index = ['Sc','Re','difl','kf','K','1/n','dp','ds','ssq','ebct','sf'])
        
        if self.optimize_flag:
            with pd.ExcelWriter(self.project_name+'_'+compound+'-'+self.carbon+'.xlsx') as writer: #'-'+repr(round(best_val_xn,2))

                model_data.to_excel(writer, 'model_fit')
                
                inf.to_excel(writer, 'influent')
                eff.to_excel(writer, 'effluent')
                data_tmp.to_excel(writer, 'parameters')
    
                ti.sleep(1)
            
        return compound, best_val_k, best_val_xn, ssqs, model_data
    #end kfit
    
    #begin dsfit
    def run_psdm_dsfit(self, compound):
        '''
        time: must be specified in minutes
        flow rate assumed to be in 'ml/min' in column properties
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
        
        k_val = self.k_data[compound]['K']
        xn = self.k_data[compound]['1/n']
        
        mw = self.comp_df[compound]['MW']            #molecular weight
        mol_vol = self.comp_df[compound]['MolarVol'] # Molar volume
        
        if self.time_type == 'days':
            dstep = 0.25 * min_per_day
        else:
            dstep = 15.        
        
        #set up bindings for nonlocal varaibles
        cinf = 1.
        cout_f = 1.
        tconv = 1.
        time_dim = 1.
        time_dim2 = 1.
        ttol = 1.
        tstep = 1.
        ds_v = 1.
                
        inf = self.data_df[self.influent][compound]
        eff = self.data_df[self.carbon][compound] 
        #convert cbo to molar values
        cbo = inf * self.mass_mul / mw #/ 1000. #(* self.mass_mult ????)
        time = (inf.index * t_mult).values
        if cbo.iloc[0] == 0.:
            cb0 = 1.
        else:
            cb0 = cbo[0]
        
        cin = cbo/cb0
        
        try:
            brk = self.k_data[compound]['brk']
        except:
            brk = np.max(inf.index.values)
               
        tortu = 1.0                             # tortuosity
        psdfr = 5.0                             # pore to surface diffusion ratio
        nd = nc - 1
        
        difl = 13.26e-5/(((vw * 100.)**1.14)*(mol_vol**0.589)) #vb
        sc = vw / (dw * difl)       #schmidt number
        
        #set film and pore diffusion
        multi_p = difl/(2*rad) # multiplier used for kf calculation
        kf_v = kf_calc(multi_p, self.re, sc, ebed, corr='Chern and Chien')

        if compound == 'Test':
            kf_v = self.k_data['Test']['kf'] #will break
            
        cout = eff * self.mass_mul / mw / cb0 
        
        dp_v = (difl/(tortu))       #porosity not used in AdDesignS appendix, removed to match
        ds_base = 1. #set up for nonlocal
        
        # @stopit.threading_timeoutable()
        def run(ds_mult):
            nonlocal cinf
            nonlocal cout_f
            nonlocal tconv
            nonlocal time_dim
            nonlocal time_dim2
            nonlocal ttol
            nonlocal tstep
            nonlocal ds_v # for output
            nonlocal ds_base
            
            aau = np.zeros(mc)
            
            #==============================================================================
            # #converts K to (umole/g)(L/umole)**(1/n), assumes units of (ug/g)(L/ug)**(1/n)
            #==============================================================================
            molar_k = k_val / mw / ((1. / mw) ** xn)  
            xni = 1./xn
            
            ds_base = epor*difl*cb0*psdfr/(1e3*rhop*molar_k*cb0**xn) 
            if self.optimize_flag:
                ds_v = ds_base * ds_mult
            else:
                ds_v = ds_mult
            
            #multiplies ds by ds_mult, passed as argument
            d = ds_v/dp_v
            
            qe = molar_k * cb0**xn 
            qte = 1. * qe
            
            dgs = (rhop * qe * (1.-ebed) * 1000.)/(ebed * cb0)
            dgp = epor * (1. - ebed)/(ebed) 
            dg = dgs + dgp
            dgt = dg
            dg1 = 1. + dgt
            dgI = 1.0/dg
            edd = dg1/dg #dgt changed from dg1
            
            ym = qe/qte
            
            eds = ds_v*dgs*tau/(rad**2)
            if eds < 1e-130:
                eds = 1e-130
            edp = dp_v*dgp*tau/(rad**2)
            
            # from orthog(n)
            beds = (eds + d*edp) * edd * br[:-1]
            bedp = edp * (1. - d) * edd * br[:-1]
            
            #depends on kf
            st = kf_v * (1. -ebed) * tau/(ebed*rad)
            stdv = st * dgt     # dgt changed from dg1
            
            #convert to dimensionless
            tconv = 60./(tau*dg1)
            tstep = dstep * tconv
            ttol = time[-1] * tconv
            time_dim = time * tconv
            
            numb = int(brk*2 + 1)
            time_dim2 = np.linspace(0., brk * tconv * t_mult,\
                                    num=numb, endpoint=True) #increase the number of sites to check for ssq
            
            #set up time based influent data
            cinf = interp1d(time_dim, cin.values, fill_value='extrapolate') 
            cout_f = interp1d(time_dim, cout.values, fill_value='extrapolate')
            #initialize storage arrays/matrices
            n = (nc+1)*mc
            y0 = np.zeros(n)

            def diffun(t, y0):
                nonlocal aau
                y0tmp = y0.reshape((nc+1,mc))
                ydot = np.zeros(y0tmp.shape)
                
                #defines the influent concentration at time t
                cinfl = cinf(t)
                fac = 1.
                
                z = ym * y0tmp[:nc,:mc] #* ym #updated ym should always be 1 for single comp.
                qte = z
                yt0 = xni * z
                
                z_c = z/qte
                z[qte>0.] = z_c[qte>0.] # should be 1 for single component.
                q0 = yt0 * xn/ym
                
                q0[np.logical_not(np.isfinite(q0))] = 0.
                z[np.logical_not(np.isfinite(z))] = 0.
                
                cpore = z * q0**xni 
                cpore[np.logical_or.reduce((qte<=0.,\
                                            yt0<=0,\
                                            xni*np.log10(q0)<-20,\
                                            cpore==np.nan
                                            ))] = 0.
                
                cpore[np.isinf(cpore)] = 1.

                cpore_tmp = cpore[nc-1]
                cpore_tmp[cpore_tmp < 0.] = 0.
                cbs = stdv*(y0tmp[nc]-cpore_tmp)
                cbs[0] = 0. 
                
                bb = fac * np.dot(bedp, cpore) + np.dot(beds,y0tmp[:nc,:])
                ww = np.dot(wr[:nd], bb)
                ydot[:nd,:] = bb
            
                ydot[nc-1][0] = (stdv*dgI*(cinfl - cpore[nc-1][0]) - ww[0]) / wr[nc-1] #iii
                ydot[nc-1][1:] = (cbs[1:]*dgI - ww[1:])/wr[nc-1]
                
                aau[1:] = (np.dot(az[1:,1:],y0tmp[-1,1:]))
                
                ydot[-1,1:] = (-dgt*(az[:,0]*cinfl + aau) - 3.* cbs)[1:]  #dgt was changed from dg1  
                ydot = ydot.reshape((nc+1)*(mc))
                return ydot
            
            try:
                y = solve_ivp(diffun,\
                                (0, ttol),\
                                y0, \
                                method=self.solver,\
                                max_step=tstep/3,\
                                )
                # defines interpolating function of predicted effluent
                cp_tmp = y.y[-1]
                cp_tmp[cp_tmp < 0.] = 0.#sets negative values to 0.
                cp_tmp[cp_tmp > np.max(cin)*3.] = np.max(cin)*3. #sets the max to 3x cb0
                cp = interp1d(y.t, cp_tmp, fill_value='extrapolate') 
            except Exception:# as e:
                t_temp = np.linspace(0, ttol, 20)
                cp_tmp = np.zeros(20)
                cp = interp1d(t_temp, cp_tmp, fill_value='extrapolate')
            return cp
        
        def run_fit(ds_mult):
            #passes a ds_multiplier, rather than ds directly
            cp = run(ds_mult)
            ssq = ((cout_f(time_dim2)-cp(time_dim2))**2).sum()
            return ssq
        
        if self.optimize_flag:
            #reuse test_range
            test_range = 10**(1-self.test_range)
            ssqs = pd.Series([run_fit(i) for i in test_range], \
                                  index=test_range)
            min_val = ssqs[ssqs==ssqs.min()]
            best_val_ds = min_val.index[0] * ds_base
            
            best_fit = run(min_val.index[0])
            min_val = min_val.values[0]
            
            with pd.ExcelWriter('ssq_'+self.carbon+'-'+compound+'.xlsx') as writer:
                ssqs.to_excel(writer, 'Sheet1')
        
        else: #assume test_range and xn_range are single values
            best_val_ds = self.test_range[0] * ds_base
            best_fit = run(best_val_ds)
            min_val = 1e2 #run_fit(best_val_k, best_val_xn)
            ssqs = pd.Series(min_val, index=[best_val_ds])
        
        itp = np.arange(0., ttol+tstep, tstep) 
        output_fit = interp1d(itp/tconv, \
                              best_fit(itp) * cb0 * mw / \
                              self.mass_mul, \
                              fill_value='extrapolate')
            
        with pd.ExcelWriter(self.project_name+'_'+compound+'-'+self.carbon+'.xlsx') as writer:
        
            model_data = pd.DataFrame(output_fit(itp/tconv), \
                                      columns = ['data'], \
                                      index = itp/tconv/t_mult)
            model_data.to_excel(writer, 'model_fit')
            
            inf.to_excel(writer, 'influent')
            eff.to_excel(writer, 'effluent')
            
            data_tmp = pd.Series([sc, self.re, difl, kf_v, self.k_data[compound]['K'],\
                                  self.k_data[compound]['1/n'], dp_v, \
                                  best_val_ds, min_val, self.ebct, self.sf], \
                                  index = ['Sc','Re','difl','kf','K','1/n','dp',\
                                           'ds','ssq','ebct','sf'])
            data_tmp.to_excel(writer, 'parameters')
            if self.optimize_flag:
                ti.sleep(1)
            
        return compound, best_val_ds, ssqs, model_data, ds_base
    #end dsfit
            
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

    def run_all_smart(self, plot=False, save_file=True, 
                      file_name='PSDM_', pm=0, num=11, des_xn=0.025):
        '''
        Smart Optimizer for K & 1/n fitting.
        Precalculates the effective fouling based on breakthrough time.
        1. Performs one scan for range of 1/n (xn_range) at estimated fouled K.
        2. Finds minimum for 1/n, then scans an additional K range (based on pm value).
        3. If further minimization is found after step 2, then an additional small
           K & 1/n polishing step is performed.
        
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

        Returns
        -------
        None.
        Best K & 1/n values are stored to self.k_data

        '''
        
        def min_fun(x0, compound, k_val, q_val, status_print=False):
            ''' x0 = [K, 1/n] '''
            kfact, xns = x0
            base_k = recalc_k(k_val, q_val, self.xn, xns)
            self.optimize_flag = False
            
            self.test_range = np.array([kfact * base_k])
            self.xn_range = np.array([xns])
            
            _, _, _, ssqs, _ = self.run_psdm_kfit(compound)
            
            if status_print:
                output = ssqs.values[0][0]
                if output > 1e-3:
                    printer = np.round(output,3)
                else:
                    printer = np.format_float_scientific(output, precision=3)
                print(np.round(ssqs.index[0],3), '\t',
                      np.round(ssqs.columns[0],3),'\t: ', 
                      printer)
            
            return ssqs.values[0][0]
        
        opt_flg = self.optimize_flag
        orig_test_range = self.test_range * 1.
        orig_xn_range = np.round(self.xn_range * 1.,5)
        self.optimize_flag = False
        
        file_name = file_name.replace(' ', '_') #removes spacing for underscore
        if file_name[-1] != '-' or file_name[-1] != '_':
            file_name += '_' #adds spacer underscore to filename
        
        for compound in self.compounds:
            print(compound, ' running')
            k_val = self.k_data[compound]['K']
            q_val = self.k_data[compound]['q']
            brk_day = self.k_data[compound]['brk'] * self.t_mult
            
            comp_fouling = self.fouling_dict[compound]
            integral, _ = quad(comp_fouling, 0, brk_day) #ignore error from quad
            k_factor = brk_day/integral
            
            #start in middle of range, xn=0.45
            best_xn = 0.45
            best_k_factor = k_factor * 1.
            x0 = [k_factor, 0.45]
            best_ssq = min_fun(x0, compound, k_val, q_val)
            
            #sets up bounds on searches
            max_xn = np.max(orig_xn_range) + 1e-4 
            min_xn = np.min(orig_xn_range) - 1e-4
            perc_pm = pm/100.
            max_k_factor = k_factor * (1 + perc_pm)
            min_k_factor = k_factor * (1 - perc_pm)
            
            #test xn range
            decreasing = True
            correct_direction = False
            sign = 1
            xn = 0.45 + sign * des_xn # initiates xn
            count = 0
            while xn <= max_xn and decreasing:
                x0 = [k_factor, xn]
                ssq = min_fun(x0, compound, k_val, q_val)
                
                if ssq < best_ssq:
                    best_xn = xn * 1
                    #best_k_factor doesn't change
                    best_ssq = ssq * 1
                    decreasing = True
                    count += 1
                    correct_direction = True
                    
                if decreasing and ssq > best_ssq:
                    decreasing = False
                
                xn = np.round(xn + sign * des_xn, 4) #increment xn
            
            if not correct_direction:
                sign = -1
                xn = 0.45 + sign * des_xn # initiates xn
                count = 0
                decreasing = True
                while xn >= min_xn and decreasing:
                    x0 = [k_factor, xn]
                    ssq = min_fun(x0, compound, k_val, q_val)
                    
                    if ssq < best_ssq:
                        best_xn = xn * 1
                        #best_k_factor doesn't change
                        best_ssq = ssq * 1
                        decreasing = True
                        count += 1
                        
                    if decreasing and ssq > best_ssq:
                        decreasing = False
                    
                    xn = np.round(xn + sign * des_xn, 4) #increment xn
            
            # increment k_range, assume best_xn is correct
            if pm > 0:
                if num==0:
                    #prevents div/zero error
                    num=1
                des_k = (max_k_factor - k_factor) / num
                
                decreasing = True
                correct_direction = False
                sign = 1
                pmk = best_k_factor + sign * des_k # initiates xn
                count = 0
                while pmk <= max_k_factor and decreasing:
                    x0 = [pmk, best_xn]
                    ssq = min_fun(x0, compound, k_val, q_val)
                    
                    if ssq < best_ssq:
                        best_k_factor = pmk * 1
                        best_ssq = ssq * 1
                        decreasing = True
                        count += 1
                        correct_direction = True
                        
                    if decreasing and ssq > best_ssq:
                        decreasing = False
                    
                    pmk += sign * des_k
                
                if not correct_direction:
                    sign = -1
                    pmk = best_k_factor + sign * des_k # initiates xn
                    count = 0
                    decreasing = True
                    while pmk >= min_k_factor and decreasing:# and pmk > 0:
                        x0 = [pmk, best_xn]
                        ssq = min_fun(x0, compound, k_val, q_val)
                        
                        if ssq < best_ssq:
                            best_k_factor = pmk * 1
                            best_ssq = ssq * 1
                            decreasing = True
                            count += 1
                            
                        if decreasing and ssq > best_ssq:
                            decreasing = False
                        
                        pmk += sign * des_k
                        
                # finer step-size polish - End step
                #search small 1/n range to see if it moved
                # print('Entering Polishing Step')
                stps = 6 #number of steps to try in 1/n range
                new_xn_r = np.linspace(best_xn-stps*des_xn,
                                        best_xn+stps*des_xn,
                                        2*stps+1)
                new_xn_range = [i for i in new_xn_r if (i <= np.max(orig_xn_range) and i >= np.min(orig_xn_range))]
                
                #currently this creates standard incriment of 0.01 or half
                #of des_k, but should this just be des_k/10.
                if des_k > 0.01:
                    des_k = 0.01
                else:
                    des_k /= 2.
                
                decreasing = True
                correct_direction = False
                sign = 1
                pmk = best_k_factor + sign * des_k # initiates xn
                
                while pmk <= max_k_factor and decreasing:
                    count = 0

                    for xn in new_xn_range:
                        x0 = [pmk, xn]
                        ssq = min_fun(x0, compound, k_val, q_val)

                        if ssq < best_ssq:
                            best_xn = xn * 1
                            best_k_factor = pmk * 1
                            best_ssq = ssq * 1
                            decreasing = True
                            count += 1
                            correct_direction = True
                            
                        if decreasing and ssq > best_ssq and count==0:
                            # it means that no 1/n value caused the ssq to decrease
                            decreasing = False
                    
                    pmk += sign * des_k
                
                if not correct_direction:
                    sign = -1
                    pmk = best_k_factor + sign * des_k # initiates xn
                    decreasing = True
                    
                    while pmk >= min_k_factor and decreasing:
                        count = 0
                        for xn in new_xn_range:
                            x0 = [pmk, xn]
                            ssq = min_fun(x0, compound, k_val, q_val)

                            if ssq < best_ssq:
                                best_xn = xn * 1
                                best_k_factor = pmk * 1
                                best_ssq = ssq * 1
                                decreasing = True
                                count += 1
                                
                            if decreasing and ssq > best_ssq and count==0:
                                decreasing = False
                        
                        pmk += sign * des_k
           
            #convert back to best_k
            best_k = best_k_factor * recalc_k(k_val, q_val, self.xn, best_xn)
            
            print(f"K: {best_k:.3f} -- 1/n: {best_xn:.3f}")
            self.k_data[compound]['K'] = best_k * 1
            self.k_data[compound]['1/n'] = best_xn * 1
            
            if plot: 
                inf = self.data_df[self.influent][compound]
                eff = self.data_df[self.carbon][compound]
                
                plt.plot(inf.index, inf.values, marker='.', ls=':', 
                         color='silver', 
                         markerfacecolor='None', label='influent') #empty circle
                plt.plot(eff.index, eff.values, marker='+', ls='None', 
                         color='grey', label='effluent')
                
                #re-run best fit... 
                self.test_range = np.array([best_k])
                self.xn_range = np.array([best_xn])
                _, _, _, _, md = self.run_psdm_kfit(compound)
                
                plt.plot(md.index, md.values, ls='-',marker='None',color='black',
                         label='model')
                
                plt.title(f"{compound}\nK: {best_k:.3f} - 1/n: {best_xn:.3f}")
                plt.xlabel('Time (days)')
                plt.ylabel('Concentration (ng/L)')
                plt.legend()
                plt.savefig(file_name+compound+'.png',dpi=300)
                plt.close()
            
        if save_file:
            '''
            need to add save file handler. 
            currently not implimented
            '''
            # writer = pd.ExcelWriter('best_fits-'+self.project_name+'.xlsx')
            with pd.ExcelWriter('best_fits-'+self.project_name+'.xlsx') as writer:
                self.k_data.to_excel(writer, 'Sheet1')
        
        #reset original values
        self.optimize_flag = opt_flg
        self.test_range = orig_test_range * 1.
        self.xn_range = orig_xn_range * 1.

#END RUN_ALL_SMART
    
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
        ds_v = 1. #### need to add in user input mass transfer ###TODO!
                
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
        
        
        
        # @stopit.threading_timeoutable()
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
            ym_vA = ym_v.reshape(ThreeDSize) #.__altshape
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
                                fill_value = 'extrapolate') #kind='cubic' 
                
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
                qte2[qte2==0.] = 1e-30 ### set mass to very small number (avoid divide by zero)
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
                            max_step=tstep/3.,\
                            )
            
            # defines interpolating function of predicted effluent
            time_index = y.t/tconv/t_mult # creates time_array and converts back to dimensions 

            cp_tmp = y.y[effluent_locator].transpose()
            cp_tmp[cp_tmp < 0.] = 0.#sets negative values to 0.
            
            cp_df = pd.DataFrame(cp_tmp, index=time_index, columns=compound_list)
            cp_df = cp_df * cb0 * mw / self.mass_mul
            # print(cp_df)
        
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
    #end multi_run
    
    
    
# =============================================================================
#     NEW BELOW, may delete
# =============================================================================
    
    # def __run_MP_helper(self, k, invN, compound, k_mult):
    #     mp.freeze_support()
    #     self.test_range = np.array([k])
    #     self.xn_range = np.array([invN])
    #     compound, best_val_k, best_val_xn, ssqs, model_data = self.run_psdm_kfit(compound)
    #     print(ssqs)
    #     return k, invN, ssqs, compound, k_mult
        
    
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
