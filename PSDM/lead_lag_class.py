# -*- coding: utf-8 -*-
"""
Lead-Lag Modeling Class for GAC and IX media

This code represents a proposed method for handling lead-lag systems for 
GAC-GAC, IX-IX, and GAC-IX or IX-GAC configurations. This interacts with both
PSDM and HSDMIX codes, and assumes that fouling within these systems behaves
like single column cinfigurations. No guarantee is made about the accuracy of 
these simulation results. 

@author: JBurkhar
"""

import sys
psdm_loc = 'path_to_PSDM' ## example 'C:/Users/JBurkhar/OneDrive/Desktop/PSDM'
sys.path.append(psdm_loc)


import warnings
warnings.simplefilter("ignore")
import os
curdir = os.getcwd()

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
# import matplotlib.pyplot as plt
# plt.close('all') # can be commented out, this just closes any previously generated graphics

import PSDM
from PSDM_functions import lpg, min_per_day #, gm_per_lb
from ixpy import hsdmix


ml_per_gal = lpg * 1e3
s_per_day = min_per_day * 60

testing=False

## set up dictionary for flow unit conversions
flow_converter = {'ml/min': 1, # expected unit
                  'mgd': 1e6 * ml_per_gal/min_per_day, ## million gallon per day
                  'gpm': ml_per_gal, ## gallon per minute
                  'gal/min': ml_per_gal, 
                  'lps': 1e3 * 60, ## liter per second
                  'l/s': 1e3 * 60,
                  'm3ps': 1e6 * 60, ## cubic meter per s,
                  'm3/s': 1e6 * 60,
                  'l/min': 1e3, #liter per minute,
                  }


### eliminate need for exporting files, used by HSDMIX
from io import BytesIO

def io_data(data_dict):
    '''
    

    Parameters
    ----------
    data_dict : dictionary of dataframes
        must be in the form:
            
            {'Cin': cin_df,
             'params': params_df,
             'ions': ions_df
             }


    Returns
    -------
    iex_inmem : BytesIO RAM file
        a RAM file that can be supplied to HSDMIX in place of saved file

    '''
    
    iex_inmem=BytesIO()
    with pd.ExcelWriter(iex_inmem, engine='xlsxwriter') as writer:
        for sheet_name in data_dict.keys():
            data_dict[sheet_name].to_excel(writer, sheet_name=sheet_name, index=False)
    iex_inmem.seek(0,0)
    return iex_inmem

def conv_iex_u(u, t, ions):
    ix_df = pd.DataFrame(index=np.round(t/s_per_day,2), ## rounding prevents indexing error in loop
                         columns=ions['name'])
    
    num_compounds = len(ions['name'])
    for i in range(num_compounds):
            ## only convert non-counterion concentrations.. may need to fix
        factor = 1 ## assumes meq
        if ions['units'][i] == 'ug':
            factor = ions['mw'][i] * 1e3 ## meq -> ug/L
        elif ions['units'][i] == 'ng':
            factor = ions['mw'][i] * 1e6 ## meq -> ng/L
        elif ions['units'][i] == 'mg':
            factor = ions['mw'][i] * 1 ## meq -> mg/L
        ### what about mgC, mgN units??? 
        ## TODO: >>
        
        ix_df[ions['name'][i]] = u[0, i, -1, :] * factor ### returned as meq... need to add smarter conversion
        
       
    ix_df[ix_df <= 0] = 1e-16 ## make negatives or 0's very small, numerical stability
    
    return ix_df
    
    
    
    

# =============================================================================
# Begin Class
# =============================================================================

class LLobj():
    '''
    Lead Lag Handling Object Class
    
    Assumes fouling, or impact of NOM is consistent between lead and lag vessel.
    
    Can handle GAC-> GAC, IX->IX, GAC->IX, IX->GAC
    
    These fuctions are proposed, and not yet tested. No gaurantee of results 
    should be inferred.
    
    NOTE: IX:  current code assumes counterions are provided as meq, and primary
               ions are provided as mass (ng)
          
    
    
    '''
    
    def __init__(self, **kw):
        self.error = False
        self.title = kw.get('title', 'Lead Lag System - Python Generated File')
        self.job_name = kw.get('job', 'LL_file')
        
        
        self.configuration = kw.get('configuration', 'gac') ### might be redundant with LL_type... 
        # 'gac', 'ix', 'gac-ix', 'ix-gac' Viable options
        
        
        ## read in chemical properties
        self.chem_properties_file = kw.get('properties_file', 'PFAS_properties.xlsx')
        self.chem_properties = kw.get('properties', 'None')
        if type(self.chem_properties) == str:
            try:
                self.chem_properties = pd.read_excel(self.chem_properties_file, index_col=0, header=0, sheet_name='Sheet1')
            except Exception as e:
                if self.chem_properties == 'None':
                    self.error = True
                    print(f'Error reading {self.chem_properties_file}: {e}')
                    self.chem_properties = pd.DataFrame(columns=['None'], index=['MW','MolarVol','Density','Solubility','VaporPress','BP'])
        
        ## read in raw data
        self.raw_data = kw.get('raw_data', 'None')
       
        if type(self.raw_data) == str:
            self.error = True
            print('GAC: Please provide raw_data as dataframe')
        elif type(self.raw_data) != pd.core.frame.DataFrame:
            self.error = True
            print('GAC: Please provide raw_data as dataframe')  ## could be set up to support dictionaries?

        ## TODO: Should there be a way to specify GAC vs IX raw data?

        ### figure out available compounds, or be supplied list to run
        self.compounds = kw.get('compounds', 'use raw_data') ## provided as a list, or use raw_data
        if not self.error and self.compounds == 'use raw_data':
            try:
                self.compounds = list(self.raw_data['influent'].columns)    ## expected to be multiindex dataframe with influent in first level
            except:
                self.error = True
                print('Error determining available compounds, check raw_data')

        ### get GAC Freundlich isotherm data
        self.k_data = kw.get('k_data', 'None')
        
        if type(self.k_data) == str and 'gac' in self.configuration.lower():
            self.error = True
            print('GAC: No K data provided, unable to run GAC simulation')
        elif type(self.k_data) != pd.core.frame.DataFrame and 'gac' in self.configuration.lower():
            self.error = True
            print('GAC: Please provide provide K_data as a dataframe')  ## could be set up to support dictionaries?
        
        
        ## system flow input??
        #### TODO -- if gac_ix or ix_gac, the flow rates must be matched!!!
        
        ### get GAC column data or create default
        ## default configuration 10ft diameter bed with 20,000 lbs F400, 0.55 gm/ml, 10 min EBCT 
        self.column_info = kw.get('column_info', pd.DataFrame(data=[0.0513, 1.6494268e6, 0.641,    5,  0.803, 0.55, 226.1, 9.0718474e6, 304.8, 1, 'influent', 'F400', 'ng', 'days', 1, 1440, 1e-3, 'ml' ],
                                                              index=['rad', 'flrt', 'epor', 'psdfr', 'rhop', 'rhof', 'L', 'wt', 'diam', 'tortu', 'influentID','effluentID','units','time', 'mass_mul','t_mult','flow_mult','flow_type'],
                                                              columns=['F400']))

        self.influentID = self.column_info.loc['influentID']
        
        self.effluentID = self.column_info.loc['effluentID']
        # self.carbon = self.column_info.columns[0] ### should be the same as effluentID, need to check
        ## need to add something to correct the carbon name in raw_data to match with that is in effluentID in column_info
        ### TODO !!!!
        
        ## collocation points, GAC
        self.gac_nr = kw.get('gac_nr', 8)
        self.gac_nz = kw.get('gac_nz', 13)
        self.gac_ne = kw.get('gac_ne', 1)

        #########
        ## parameters needed for GAC Model
        self.water_type = kw.get('water_type', 'Organic Free')       ### no fouling
        self.chem_type = kw.get('chem_type', 'halogenated alkenes')  ### no fouling
        
        if not self.error:  ## should ensure that raw_data is a dataframe
            self.ix_raw_data = self.raw_data[self.influentID]
            self.ix_raw_data['Time'] = self.ix_raw_data.index ### will the IX code like if this isn't the 1st column...
        
        # =============================================================================
        #         ##### IX/RESIN parameters
        # =============================================================================
        self.ix_data = {'Cin': [], 'params': [], 'ions': []}
        self.ix_file = kw.get('ix_filename', 'None')  ### do we set this up as a default file, could reduce some of the other defaults below ##TODO:??
                ## collocation points, IX
        self.ix_nr = kw.get('ix_nr', 8)
        self.ix_nz = kw.get('ix_nz', 12)
        self.ix_ne = kw.get('ix_ne', 1) 
        
        if self.ix_file == 'None':
            ## default selectivies for RESIN: ----- ### TODO !!!!
            self.counterions = kw.get('counter ions', [['CHLORIDE',    35.45,     1, 1, 'meq'],
                                                       ['SULFATE',     96.06, 0.028, 2, 'meq'],
                                                       ['BICARBONATE', 61.02, 0.392, 1, 'meq'], ## 61.02, total alkalnity 50.0
                                                       ['NITRATE',        62,    13, 1, 'meq']])
            
            self.counterion_names = [i[0] for i in self.counterions]

            self.ions = list(self.counterions) ## not sure if list is needed, copy counterions to base ions 
            self.avail_ions = [self.ions[i][0] for i in range(len(self.ions))]
            
            ## params
            self.ix_params = kw.get('ix_params', pd.DataFrame(columns=['name','value','units'], data=[['Qf',1000.,'meq/L'],
                                                                                                      ['EBED',0.35, ''],
                                                                                                      ['L', 350/(np.pi/4 * 10**2), 'ft'], ## 350 ft^3 of media
                                                                                                      ['flrt', 4.9554481536e6, 'ml/min'],  ### 2 min EBCT
                                                                                                      ['rb', 0.0345, 'cm'],
                                                                                                      ['kL', 4.7e-3, 'cm/s'],
                                                                                                      ['Ds', 5e-8, 'cm2/s'],  # for counterions, for PFAS ~5e-10
                                                                                                      ['nr', self.ix_nr, ''],
                                                                                                      ['nz', self.ix_nz, ''],
                                                                                                      ['time', 1, 'day'],
                                                                                                      ['diam', 10, 'ft'],
                                                                                                      ]))
            
            ##  selectivity ### similar to k_data? ###TODO
            self.selectivities = kw.get('selectivity', {'PFBA': 127,
                                                        'PFPeA': 513,
                                                        'PFHxA': 1426,
                                                        'GenX': 1426, ##est
                                                        'PFHpA': 3311,
                                                        'PFOA': 12702,
                                                        'PFBS': 22596,
                                                        'PFNA': 22863,
                                                        'PFPeS': 69000,  ##est
                                                        'PFHxS': 116935,
                                                        'PFOS': 2e6})
            
            ### fouling adjustment for IX as [lower Kxc, upper Kxc, factor below upper Kxc, factor above upper Kxc]
            self.adjustment_factor = kw.get('selectivity adjustment', [100, 3500, 1, 1] ) ## defaults 0, 0
            
            for comp in self.selectivities.keys(): ## doing this differently later #### TODO, need to look at raw_data??
                if comp in self.compounds:
                    self.add_resin_selectivity(comp, self.selectivities[comp])
            
            ### storing things into the appropriate dictionary
            self.ix_data['params'] = self.ix_params
            self.ix_data['Cin'] = self.ix_raw_data
            self.ix_data['ions'] = self.ions
            
        else:
            ## if file name
            try:
                self.ix_data['params'] = pd.read_excel(self.ix_file, sheet_name='params', header=0)
                nr_idx = self.ix_data['params'][self.ix_data['params']['name'] == 'nr'].index.values[0]
                nz_idx = self.ix_data['params'][self.ix_data['params']['name'] == 'nz'].index.values[0]
                
                ### overrides values provided in input file --- may need to change
                self.ix_data['params'].loc[nr_idx, 'value'] = self.ix_nr
                self.ix_data['params'].loc[nz_idx, 'value'] = self.ix_nz
                
                ## NE?
                
            except:
                self.error = True
                print(f"Error reading in params from {self.ix_file}")
            try:
                self.ix_data['Cin'] = pd.read_excel(self.ix_file, sheet_name='Cin', header=0)
                
                
                time_units = (self.ix_data['params'][self.ix_data['params']['name']=='time'])['units'].values[0]
                if time_units == 'hr': ## correct time if needed
                    time_idx = self.ix_data['params'][self.ix_data['params']['name']=='time'].index.values[0]
                    self.ix_data['params'].loc[time_idx, 'units'] = 'day'
                    self.ix_data['Cin']['Time'] /= 24
                
                self.ix_raw_data = self.ix_data['Cin'].copy()    
                print('IX: Found influent concentrations')
                self.error = False ### need better handling of errors ## TODO:  error, GAC/IX split?
            except:
                self.error = True
                print(f'Error reading in Cin from {self.ix_file}')
            
            try:
                self.ix_data['ions'] = pd.read_excel(self.ix_file, sheet_name='ions', header=0)
            except:
                self.error = True
                print(f"Error reading in ions from {self.ix_file}")
            pass
        
            self.avail_ions = [self.ix_data['ions']['name'][i] for i in self.ix_data['ions'].index]
            self.counterion_names = ['CHLORIDE', 'SULFATE', 'BICARBONATE', 'NITRATE']  ### TODO: Should this be hardcoded? can we make dynamic
        
        # =============================================================================
        #         ## testing
        # =============================================================================
        if testing: 
            print(self.raw_data)
            print(self.k_data)
            print(self.compounds)
            print(self.ix_raw_data)  
        
        
        
        #### set up default results structure
        ### currently not used: May be different between IX/GAC. ##TODO: CHECK
        self.results = {'mid': {},
                        'mid_results': {},
                        'effluent': {},
                        'effluent_results': {}}
        
    def add_resin_selectivity(self, comp, selectivity):
        if comp not in self.avail_ions:  ## add if not in list
            self.avail_ions.append(comp)
            
            Kxc = self.selectivities[comp]
            if Kxc >= self.adjustment_factor[0] and Kxc < self.adjustment_factor[1]:
                Kxc *= self.adjustment_factor[2]  ## earlier eluting
            elif Kxc > self.adjustment_factor[1]:
                Kxc *= self.adjustment_factor[3]  ## later eluting
            
            self.ions.append([comp, self.chem_properties[comp]['MW'], Kxc, 1, 'ng'])   #assumes always ng, might need to change ### TODO
        else: ## update
            for i in range(len(self.ions)):
                if self.ions[i][0] == comp:
                    Kxc = self.selectivities[comp]
                    if Kxc >= self.adjustment_factor[0] and Kxc < self.adjustment_factor[1]:
                        Kxc *= self.adjustment_factor[2]  ## earlier eluting
                    elif Kxc > self.adjustment_factor[1]:
                        Kxc *= self.adjustment_factor[3]  ## later eluting
                    
                    self.ions[i][2] = Kxc * 1
        
    def update_ion_unit(self, comp, unit):
        
        pass
    
    ## primary function, calls sub runners
    def run_LL(self, swap_days, LL_type='gac', loops=5, **kw):
        compounds = kw.get('compounds', self.compounds)
        
        if self.error:
            print('There is an error in an input, please retry')
            data_store = {'error': 'There was an error'}
        else:
            if LL_type.lower() == 'gac' or LL_type.lower() == 'gac-gac':
                data_store = self.run_LL_gac(swap_days, loops=loops, compounds=compounds)
            elif LL_type.lower() == 'ix' or LL_type.lower() == 'ix-ix' or\
                 LL_type.lower() == 'iex' or LL_type.lower() == 'iex-iex':
                data_store = self.run_LL_ix(swap_days, loops=loops, compounds=compounds)
                
            elif LL_type.lower() == 'ix-gac' or LL_type.lower() == 'iex-gac':
                data_store = self.run_LL_ix_gac(swap_days, compounds=compounds, **kw)
            elif LL_type.lower() == 'gac-ix' or LL_type.lower() == 'iex-gac':
                data_store = self.run_LL_gac_ix(swap_days, compounds=compounds, **kw)
            else:
                print('You have specified a configuration that cannot be run:\nPlease use- gac, ix, gac-ix, or ix-gac')
                data_store = []
        
        return data_store
    
    ## individual runners
    ## GAC-GAC
    def run_LL_gac(self, swap_days, loops=5, **kw):
        '''
        GAC Lead - GAC Lag

        Parameters
        ----------
        swap_days : TYPE
            DESCRIPTION.
        
        compounds : list
            list of compounds to model, default is all compounds provided to class
        
        loops : int
            Sets number of loops to run, default is 5

        Returns
        -------
        None.

        '''
        compounds = kw.get('compounds', self.compounds)

        run_duration = int(swap_days * 2.1) ## make run duration slightly longer than double the swap days
        
        ### TODO: Can/Should this be moved to a special function? Can this be shared?
        raw_data_orig = self.raw_data.copy() ## make a copy of the original raw data 
        
        ## Create new influent structure
        raw_data_new = pd.DataFrame(index=np.arange(0, run_duration, 0.25), ### create a timeseries complicated enough to accept complex curves
                                    columns=raw_data_orig.columns) ## can reuse the column headers
        
        ## make everything initially the 'mean' of raw_data_orig
        for idx in raw_data_orig.columns:
            raw_data_new[idx] = raw_data_orig[idx].mean()
            
        ### create data_store object that will be returned
        data_store = {'mid':      {},
                      'mid results': {},
                      'effluent': {},
                      'effluent results':  {}
                      }
        
        for comp in compounds: ## might be able to make this loop simpler
            column = PSDM.PSDM(self.column_info,
                               self.chem_properties,
                               raw_data_new,
                               nr=self.gac_nr,
                               nz=self.gac_nz,
                               ne=self.gac_ne,
                               chem_type=self.chem_type,
                               water_type=self.water_type,
                               k_data=self.k_data,
                               optimize=False)
            
            ## set K's and 1/n's
            column.xn_range = np.array([self.k_data[comp]['1/n']])
            column.test_range = np.array([self.k_data[comp]['K']])
           
            _, _, _, _, results = column.run_psdm_kfit(comp) ## run the model
            
            for i in range(loops+1):
                raw_df2 = raw_data_new.copy()
                
                if i == 0:
                    raw_df2[self.influentID, comp].loc[:swap_days] = results['data'][:swap_days].values
                else:
                    # creates the first half of time as the effluent from the second column, replicating lead/lag
                    raw_df2[self.influentID, comp].loc[:swap_days] = results['data'][swap_days:2*swap_days].values
                
                column = PSDM.PSDM(self.column_info,
                                   self.chem_properties,
                                   raw_df2, ## use new df
                                   nr=self.gac_nr,
                                   nz=self.gac_nz,
                                   ne=self.gac_ne,
                                   chem_type=self.chem_type,
                                   water_type=self.water_type,
                                   k_data=self.k_data,
                                   optimize=False)
                
                ## set K's and 1/n's
                column.xn_range = np.array([self.k_data[comp]['1/n']])
                column.test_range = np.array([self.k_data[comp]['K']])
                
                _, _, _, _, results = column.run_psdm_kfit(comp) ## run the model
        
                if i == loops:
                    ## save final results
                    
                    ## should data_store just be stored to self.? ### TODO: XXX
                    data_store['mid'][comp] = results.loc[swap_days*2].values[0]
                    data_store['effluent results'][comp] = results.loc[0:swap_days]  ## should this just be saved as a dataframe, might be cleaner
                    data_store['effluent'][comp] = results.loc[swap_days].values[0]
                    results.index = results.index - swap_days
                    data_store['mid results'][comp] = results.loc[0:swap_days] ## should this just be saved as a dataframe, might be cleaner
                    
                    ############# NOT YET IMPLIMENTED ### NEEDED????
                    ### if we want mid point from lag vessel, need to revisit
                    # mid_point = int(np.ceil(nz/2))
                    # mid_lag_f = interp1d(column.yt, column.ydot[-mid_point],
                    #                    fill_value='extrapolate')
                    
                    # stored_data = float(mid_lag_f(2*swap_days))
        
        return data_store
    
    ## IX-IX
    def run_LL_ix(self, swap_days, loops=5, **kw):
        '''
        IX Lead - IX Lag

        Parameters
        ----------
        swap_days : TYPE
            DESCRIPTION.
        loops : TYPE, optional
            DESCRIPTION. The default is 5.
        **kw : TYPE
            DESCRIPTION.

        Returns
        -------
        data_store : TYPE
            DESCRIPTION.

        '''
        compounds = kw.get('compounds', self.compounds) ### currently running everything, could filter for compounds 
        
        data_store = {'mid':      {},
                      'mid results': {},
                      'effluent': {},
                      'effluent results':  {}
                      }
        
        run_duration = int(swap_days * 2.1) ## slightly longer than double swap_days
        tvals = np.arange(0, run_duration, 0.25)
        
        raw_data_orig = self.ix_raw_data.copy()
        
        ## opportunity to make it faster by only running single species + counterions??? ##TODO:??
        ## because most other species adsorb at low solid loading, this may be an acceptable assumption
        
        raw_data_new = pd.DataFrame(index=tvals, ### create a timeseries complicated enough to accept complex curves
                                    columns=raw_data_orig.columns) ## can reuse the column headers
        
        ## make everything initially the 'mean' of raw_data_orig
        ## variable influent doesn't make as much sense for lead-lag for optimization
        for idx in raw_data_orig.columns:
            if idx != 'Time': ## don't need to add mean to Time
                raw_data_new[idx] = raw_data_orig[idx].mean()
            else:
                raw_data_new['Time'] = tvals * 1
        
        self.ix_data['Cin'] = raw_data_new.copy()
        
        if testing:
            print(self.ix_data['Cin'])
        
        ## Run IX code
        IXcol = hsdmix.HSDMIX(io_data(self.ix_data)) ## returns simulated file in memory that can by used by IX code.
        t, u = IXcol.solve(t_eval=tvals, const_Cin=False)
        
        converted_df = conv_iex_u(u, t, self.ix_data['ions']) ## only work with effluent concentrations
        
        if testing:
            converted_df.plot.line()
            print(raw_data_new)
            print(converted_df)
            
        ### need to actually add loops section --- 
        ## TODO::: here we are
        
        for i in range(loops + 1):
            data_df2 = raw_data_new.copy()
            for comp in converted_df:
                if i == 0: ## first pass
                    data_df2[comp].loc[:swap_days] = converted_df[comp][:swap_days].values
                else:
                    if testing:
                        print(converted_df[comp][swap_days:2*swap_days+1])
                    data_df2[comp].loc[:swap_days] = converted_df[comp][swap_days:2*swap_days].values
        
            self.ix_data['Cin'] = data_df2.copy()

            ## run IX code
            IXcol = hsdmix.HSDMIX(io_data(self.ix_data)) ## returns simulated file in memory that can by used by IX code.
            t, u = IXcol.solve(t_eval=tvals, const_Cin=False)
            
            converted_df = conv_iex_u(u, t, self.ix_data['ions']) ## only work with effluent concentrations
            
            if testing:
                print(data_df2)
                converted_df.plot.line()
            
            if i == loops:
                for comp in converted_df.columns:
                    data_store['mid'][comp] = converted_df[comp][swap_days*2]
                    data_store['effluent results'][comp] = converted_df[comp].loc[0:swap_days]  ## should this just be saved as a dataframe, might be cleaner
                    data_store['effluent'][comp] = converted_df[comp][swap_days]
                
                converted_df.index = converted_df.index - swap_days
                for comp in converted_df.columns:
                    data_store['mid results'][comp] = converted_df[comp].loc[0:swap_days] ## should this just be saved as a dataframe, might be cleaner
        
        return data_store
    
    def run_LL_gac_ix(self, swap_days, **kw):
        '''
        GAC Lead - IX Lag
        
        Assumes GAC has no impact on counterions, and they pass through to IX 
        at influent conentrations
        
        Influent to GAC (raw_data) is assumed to be the influent to system
        self.ix_data['Cin'] is IGNORED
        
        swap_days = duration of simulation for this

        Parameters
        ----------
        duration : TYPE
            DESCRIPTION.
        **kw : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        
        compounds = kw.get('compounds', self.compounds)
        shared_flow = kw.get('shared_flowrate', 'use gac')
        
        ### adjust for sharing flowrate through both systems
        flow_index = self.ix_data['params'][self.ix_data['params']['name'] == 'flrt'].index[0]
        if type(shared_flow) != str:
            if type(shared_flow) != list:
                #### assume [flow, units]
                if len(shared_flow) == 2:
                    flow, units = shared_flow
                    self.ix_data['params'].loc[flow_index]['value'] = flow
                    self.ix_data['params'].loc[flow_index]['units'] = units
                else:
                    print('IF shared_flow is supplied as list, please supply as [value, units]')
            elif type(shared_flow) != float:
                #### assume flow in ml/min ?? is this the best units
                
                    self.ix_data['params'].loc[flow_index]['value'] = shared_flow
                    self.ix_data['params'].loc[flow_index]['units'] = "ml/min"
                
        elif shared_flow == 'use gac':  ## reset IX flow rate to match GAC system
            self.ix_data['params'].loc[flow_index]['value'] = self.column_info.loc['flrt']
            self.ix_data['params'].loc[flow_index]['units'] = 'ml/min'
        elif shared_flow == 'use ix': ## reset GAC flow rate to match IX
            flow_index = self.ix_data['params'][self.ix_data['params']['name'] == 'flrt'].index[0]
            
            _, flow, unit = self.ix_data['params'].loc[flow_index]
            
            self.column_info.loc['flrt'] = flow * flow_converter[unit]
            ## TODO: TEST --- reset flow rate to self.column_info.loc['flrt'] to match (needs to be in ml/min)
        
        self._check_both() ## set up what is available in GAC/IX data
        
        ### consolidate self.shared compounds and compounds (only consider what is shared)
        use_compounds = list(set(compounds + self.shared_compounds))
        ### set up input tuple
        midx_tuple = [(i, j) for i in [self.influentID, self.effluentID] for j in use_compounds]
        midx = pd.MultiIndex.from_tuples(midx_tuple)
                

        ### TODO: Can/Should this be moved to a special function? Can this be shared?
        raw_data_orig = self.raw_data.copy() ## make a copy of the original raw data 
        
        ## Create new influent structure
        tvals = np.arange(0, swap_days+0.25, 0.25)
        tvals2 = tvals[:-1] ## some data includes all but last point
        raw_data_new = pd.DataFrame(0, index=tvals, ### create a timeseries complicated enough to accept complex curves
                                    columns=midx) ## can reuse the column headers
    
        if swap_days <= np.max(self.raw_data.index):
            for comp in raw_data_new.columns:
                f_comp = interp1d(raw_data_orig.index, raw_data_orig[comp], fill_value='extrapolate')
                raw_data_new[comp] = f_comp(tvals)
        else: ## use mean instead
            for comp in raw_data_new.columns:
                raw_data_new[comp] = raw_data_orig[comp].mean()
        
        # raw_data_new[np.isnan(raw_data_new)] = 0 ## replace nan's with 0 if needed 
        
        ## create data store object dictionary
        data_store = {'mid': {},
                      'mid results': raw_data_new[self.effluentID][:swap_days-0.25] * 0.,
                      'effluent': {},
                      'effluent results':  raw_data_new[self.effluentID][:swap_days-0.25] * 0. }
        
        # print(data_store)
        
        for comp in compounds: ## might be able to make this loop simpler
            column = PSDM.PSDM(self.column_info,
                               self.chem_properties,
                               raw_data_new,
                               nr=self.gac_nr,
                               nz=self.gac_nz,
                               ne=self.gac_ne,
                               chem_type=self.chem_type,
                               water_type=self.water_type,
                               k_data=self.k_data,
                               optimize=False)
            
            ## set K's and 1/n's
            column.xn_range = np.array([self.k_data[comp]['1/n']])
            column.test_range = np.array([self.k_data[comp]['K']])
           
            _, _, _, _, results = column.run_psdm_kfit(comp) ## run the model
        
            data_store['mid results'][comp] = results.values * 1
            data_store['mid'][comp] = results.values[-1][0]
        
        ix_data = pd.DataFrame(0, index=tvals2,
                               columns=['Time']+self.counterion_names+use_compounds)
        
        ## read in lead vessel effluent concentrations as influent to lag vessel
        for comp in use_compounds:
            ix_data[comp] = data_store['mid results'][comp].values
            # ix_data[comp] = raw_data_new[self.influentID, comp].mean()
        for counterion in self.counterion_names:
            ix_data[counterion] = self.ix_data['Cin'][counterion].mean() ## can only do mean??? maybe not?
        
        ix_data[ix_data <= 0] = 1e-16 ## reset negative/zero to very small number
        ix_data['Time'] = tvals2 ## set time
        
        self.ix_data['Cin'] = ix_data.copy()
            
        if testing:
            print(data_store['mid results'].index)
            print(self.ix_data)
            
        IXcol = hsdmix.HSDMIX(io_data(self.ix_data)) ## returns simulated file in memory that can by used by IX code.
        t, u = IXcol.solve(t_eval=tvals2, const_Cin=False) #t_eval=tvals2
        
        converted_df = conv_iex_u(u, t, self.ix_data['ions']) ## only work with effluent concentrations
        for comp in converted_df.columns:
            data_store['effluent'][comp] = converted_df[comp].values[-1]
        
        data_store['effluent results'] = converted_df.copy() ## add interpolation first
    
        return data_store
        
    def run_LL_ix_gac(self, swap_days, **kw):
        '''
        IX Lead - GAC Lag
        
        
        Uses influent conentrations from IX column

        Parameters
        ----------
        swap_days : TYPE
            DESCRIPTION.
        **kw : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        
        
        compounds = kw.get('compounds', self.compounds)
        shared_flow = kw.get('shared_flowrate', 'use ix')
        
        ### adjust for sharing flowrate through both systems
        flow_index = self.ix_data['params'][self.ix_data['params']['name'] == 'flrt'].index[0]
        if type(shared_flow) != str:
            if type(shared_flow) != list:
                #### assume [flow, units]
                if len(shared_flow) == 2:
                    flow, units = shared_flow
                    self.ix_data['params'].loc[flow_index]['value'] = flow
                    self.ix_data['params'].loc[flow_index]['units'] = units
                else:
                    print('IF shared_flow is supplied as list, please supply as [value, units]')
            elif type(shared_flow) != float:
                #### assume flow in ml/min ?? is this the best units
                
                    self.ix_data['params'].loc[flow_index]['value'] = shared_flow
                    self.ix_data['params'].loc[flow_index]['units'] = "ml/min"
                
        elif shared_flow == 'use gac':  ## reset IX flow rate to match GAC system
            self.ix_data['params'].loc[flow_index]['value'] = self.column_info.loc['flrt']
            self.ix_data['params'].loc[flow_index]['units'] = 'ml/min'
        elif shared_flow == 'use ix': ## reset GAC flow rate to match IX
            flow_index = self.ix_data['params'][self.ix_data['params']['name'] == 'flrt'].index[0]
            
            _, flow, unit = self.ix_data['params'].loc[flow_index]
            
            self.column_info.loc['flrt'] = flow * flow_converter[unit]
            ## TODO: TEST --- reset flow rate to self.column_info.loc['flrt'] to match (needs to be in ml/min)
        
        self._check_both() ## set up what is available in GAC/IX data


        
                


        ### consolidate self.shared compounds and compounds (only consider what is shared)
        use_compounds = list(set(compounds + self.shared_compounds))
        ### set up input tuple
        midx_tuple = [(i, j) for i in [self.influentID, self.effluentID] for j in use_compounds]
        midx = pd.MultiIndex.from_tuples(midx_tuple)
                

        ### TODO: Can/Should this be moved to a special function? Can this be shared?
        raw_data_orig = self.raw_data.copy() ## make a copy of the original raw data 
        
        ## Create new influent structure
        tvals = np.arange(0, swap_days+0.25, 0.25)
        tvals2 = tvals[:-1] ## some data includes all but last point
        raw_data_new = pd.DataFrame(0, index=tvals, ### create a timeseries complicated enough to accept complex curves
                                    columns=midx) ## can reuse the column headers
        
        ### set up influent concentrations
        ix_data = pd.DataFrame(0, index=tvals,
                               columns=['Time']+self.counterion_names+use_compounds)
        
        ## read in lead vessel effluent concentrations as influent to lag vessel
        if swap_days > np.max(self.ix_data['Cin']['Time']):
            ### currently assumes days
            ## uses mean values if swap_days exceeds ix_data influent
            for comp in self.counterion_names+use_compounds:
                ix_data[comp] = self.ix_data['Cin'][comp].mean()
        else:
            for comp in self.counterion_names+use_compounds:
                f_comp = interp1d(self.ix_data['Cin']['Time'].values, 
                                  self.ix_data['Cin'][comp].values,
                                  fill_value='extrapolate')
                ix_data[comp] = f_comp(tvals)
                
        ix_data[ix_data <= 0] = 1e-16 ## reset negative/zero to very small number
        ix_data['Time'] = tvals ## set time
        
        self.ix_data['Cin'] = ix_data * 1
        
        data_store = {'mid': {},
                      'mid results': ix_data.copy,
                      'effluent': {},
                      'effluent results':  raw_data_new[self.effluentID][:swap_days-0.25] * 0. }

        ### RUN IX LEAD
        IXcol = hsdmix.HSDMIX(io_data(self.ix_data)) ## returns simulated file in memory that can by used by IX code.
        t, u = IXcol.solve(t_eval=tvals, const_Cin=False) #t_eval=tvals2
        
        converted_df = conv_iex_u(u, t, self.ix_data['ions']) ## only work with effluent concentrations
        for comp in converted_df.columns:
            data_store['mid'][comp] = converted_df[comp].values[-1]
        
        data_store['mid results'] = converted_df.copy() ## add interpolation first

        thresh = 1e-8 ## to handle errors sometimes caused by low influent conditions caused by IX effluent (run_psdm_kfit throws errors)
        converted_df[converted_df <= thresh ] = thresh ## make negative/zero very small 
        ## SET UP GAC LAG
        for comp in use_compounds:
            raw_data_new[self.influentID, comp] = converted_df[comp]
        
        if testing:
            print(raw_data_new.plot.line())
            print(self.k_data)
        
        ## RUN GAC LAG
        for comp in use_compounds:
            column = PSDM.PSDM(self.column_info,
                               self.chem_properties,
                               raw_data_new,
                               nr=self.gac_nr,
                               nz=self.gac_nz,
                               ne=self.gac_ne,
                               chem_type=self.chem_type,
                               water_type=self.water_type,
                               k_data=self.k_data,
                               optimize=False)
            
            ## set K's and 1/n's
            column.xn_range = np.array([self.k_data[comp]['1/n']])
            column.test_range = np.array([self.k_data[comp]['K']])
           
            _, _, _, _, results = column.run_psdm_kfit(comp) ## run the model

            data_store['effluent'][comp] = results.values[-1][0]
            data_store['effluent results'][comp] = results.values

        return data_store
    
    def _check_both(self):
        
        self.error_both = False ## create error holder if something is missing
        ### may not need to do this...
        ### assumes defaults are provided
        
        
        
        self.shared_compounds = [i for i in self.raw_data[self.influentID].columns if i in self.ix_data['Cin'].columns]
        
        
        
        
    







