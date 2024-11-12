# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 15:45:43 2022

@author: JBurkhar
"""

# import lead_lag_class
from lead_lag_class import LLobj

import sys
psdm_loc = 'path_to_PSDM' ## example 'C:/Users/JBurkhar/OneDrive/Desktop/PSDM'
sys.path.append(psdm_loc)


import warnings
warnings.simplefilter("ignore")
import os
srt_dir = os.getcwd()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.close('all') # can be commented out, this just closes any previously generated graphics


#Comment next line in to CHANGE TO WHERE THE PSDM IS STORED ON YOUR MACHINE
# os.chdir('../')
import PSDM
os.chdir(srt_dir)

## predefine compounds included with the input files
included_compounds = ['PFBA', 'PFPeA', 'PFHxA', 'PFHpA', 'PFOA', 'PFBS', 'PFHxS', 'PFOS']

dpi = 150 ## plotting option

test_gac_gac = True 
test_ix_ix = True 
test_gac_ix = True 
test_ix_gac = True


swap_days = 130

### Generic Testing
####   TEST GAC
if test_gac_gac:
    print('################# Testing GAC-GAC ##############################')
    
    
    fn = 'test_LL/test_gac.xlsx'

    # =============================================================================
    # READ IN Chemical Property Information
    # may need to change location of file, depending on where it is stored locally
    # properties are included in 'fn' from above, but can be separate file as well
    # =============================================================================
    chem_data = PSDM.process_input_data(fn, sheet_name='Properties') 
    k_data = pd.read_excel(fn, sheet_name='Kdata',index_col=0) # K & 1/n data

    # =============================================================================
    # Read in Input File
    # Can change location of data_sheet= and column_sheet= to match file naming
    # =============================================================================
    raw_data, column_info,\
    compounds, carbons, = PSDM.process_input_file(fn,\
                                                   data_sheet='data',\
                                                   column_sheet='columnSpecs'
                                                  )
    
    LL = LLobj(raw_data=raw_data,
                   properties=chem_data,
                   column_info=column_info[carbons[0]],
                   k_data=k_data,
                   )
    
    
    if False:
        print(LL.chem_properties)
        print(LL.counterions)
        print(LL.column_info)
        print(LL.ix_params)
        print(LL.avail_ions)
        print(LL.ions)

   
    # swap_days = 300
    plt.figure(dpi=dpi)
    if False:
        # =============================================================================
        # run single column - to compare
        # =============================================================================
        for comp in compounds:
            # SETTING UP PSDM simulation column
            column = PSDM.PSDM(column_info[carbons[0]], 
                            chem_data, 
                            raw_data,
                            nr=8,
                            nz=12, 
                            k_data=k_data,\
                            optimize=False
                            )
            
            column.test_range = np.array([k_data[comp]['K']])
            column.xn_range = np.array([k_data[comp]['1/n']])
            
            # RUNNING the simulation
            #only results is used in this example
            _compound, _k, _xn, _ssqs, results = column.run_psdm_kfit(comp)
            
            
            plt.plot(results.index, 
                    results.values,
                    ls=':',
                    label=f"{comp} - Single")
        
        
    # =============================================================================
    # run lead-lag for gac    
    # =============================================================================
    LL_gac = LLobj(raw_data=raw_data,
                   properties=chem_data,
                   column_info=column_info[carbons[0]],
                   k_data=k_data,
                   water_type='Rhine',
                   chem_type='PFAS',
                   )
    
    results_dict = LL_gac.run_LL(swap_days, LL_type='gac')
    for comp in included_compounds:
        effl = results_dict['effluent results'][comp]
        p = plt.plot(effl.index, effl.values, ls='-', label=f"{comp}")
        
        if True:
            mid = results_dict['mid results'][comp]
            plt.plot(mid.index, mid.values, ls='--', alpha=0.5, color=p[0].get_color()) ## shows mid-point as dashed line
        
        
        
    
    plt.legend(loc='upper left', ncols=2, )
    plt.title(f"GAC-GAC (Swap Day: {swap_days})", fontsize=20)
    plt.ylabel('Concentration (ng/L)', fontsize=20)
    plt.xlabel('Time (days)', fontsize=20)
    plt.xlim((0, swap_days))

### Test IX
if test_ix_ix: 
    print('################# Testing IX-IX ##############################')
    plt.figure(dpi=dpi)

    ix_fn = 'test_LL/lag_iex.xlsx'
    
    LL_ix = LLobj(ix_filename=ix_fn, ix_nr=7, ix_nz=13)
    
    # swap_days = 300
    results = LL_ix.run_LL(2*swap_days, LL_type='ix')
    # print(results)
    
    for comp in included_compounds:
        temp_res = results['effluent results'][comp]
        if comp not in ['CHLORIDE','SULFATE','BICARBONATE','NITRATE']:
            p = plt.plot(temp_res.index,
                        temp_res.values, ## file had mg/L specified for PFAS
                        label=comp)
            
            plt.plot(results['mid results'][comp].index,
                     results['mid results'][comp].values,
                     alpha=0.5,
                     #f"GAC-Lead: {comp}",
                     ls=':', 
                     color = p[0].get_color())
            
    plt.legend(loc='upper left', ncols=2, title='Lag: Solid, Lead: Dotted')
    plt.title(f"IX-IX (Swap Day: {2*swap_days})", fontsize=20)
    plt.ylabel('Concentration (ng/L)', fontsize=20)
    plt.xlabel('Time (days)', fontsize=20)
    
    if False: ## turn on if you want to see major anions
        plt.figure()
        for comp in results['effluent'].keys():
            temp_res = results['effluent results'][comp]
            if comp in ['CHLORIDE','SULFATE','BICARBONATE','NITRATE']:
                mw = LL_ix.ix_data['ions'][LL_ix.ix_data['ions']['name'] == comp]['mw'].values[0]
                val = LL_ix.ix_data['ions'][LL_ix.ix_data['ions']['name'] == comp]['valence'].values[0]
                
                plt.plot(temp_res.index,
                        temp_res.values * mw / val,
                        label=comp)
        
        plt.legend(loc='upper left',)
        plt.xlim((0, 4))
        plt.title(f"IX-IX (Swap Days: {2*swap_days})", fontsize=20)
        plt.ylabel('Concentration (mg/L)', fontsize=20)
        plt.xlabel('Time (days)', fontsize=20)
    
    
#### TEST heterogeneous 
if test_gac_ix or test_ix_gac:
    

    ix_fn = 'test_LL/lag_iex.xlsx'
    gac_fn = 'test_LL/test_gac.xlsx'
    
    
    chem_data = PSDM.process_input_data(gac_fn, sheet_name='Properties') 
    k_data = pd.read_excel(gac_fn, sheet_name='Kdata',index_col=0) # K & 1/n data

    # =============================================================================
    # Read in Input File
    # Can change location of data_sheet= and column_sheet= to match file naming
    # =============================================================================
    raw_data, column_info,\
    compounds, carbons, = PSDM.process_input_file(gac_fn,\
                                                   data_sheet='data',\
                                                   column_sheet='columnSpecs'
                                                  )
    
    # swap_days = 300

    if test_gac_ix:
        print('################# Testing GAC-IX ##############################')

        LL_gac_ix = LLobj(raw_data=raw_data,
                            properties=chem_data,
                            column_info=column_info[carbons[0]],
                            k_data=k_data,
                            ix_filename=ix_fn,
                            water_type='Rhine',
                            chem_type='PFAS',
                        )
        
        results = LL_gac_ix.run_LL(2*swap_days, LL_type='gac-ix')
    
        plt.figure(dpi=dpi)
        for comp in included_compounds:
            try:
                p = plt.plot(results['mid results'].index,
                        results['mid results'][comp].values,
                        alpha=0.5,
                        #f"GAC-Lead: {comp}",
                        ls=':')
                color = p[0].get_color()
                plt.plot(results['effluent results'].index,
                        results['effluent results'][comp].values,
                        label=comp, # label=f"IX-Lag: {comp}",
                        color=color)
            except:

                pass
        
        plt.legend(loc='upper left', ncols=2, title="Solid: IX-Lag, Dotted: GAC-Lead")

        plt.title(f'GAC Lead-IX Lag (Swap Day: {2*swap_days})',  fontsize=20)
        plt.ylabel('Concentration (ng/L)', fontsize=20)
        plt.xlabel('Time (days)', fontsize=20)
        plt.tight_layout() 


    
    
    
    ### IX-GAC 
    if test_ix_gac:
        print('################# Testing IX-GAC ##############################')

        LL_ix_gac = LLobj(raw_data=raw_data,
                            properties=chem_data,
                            column_info=column_info[carbons[0]],
                            k_data=k_data,
                            ix_filename=ix_fn,
                            water_type='Rhine',
                            chem_type='PFAS',
                        )
        results2 = LL_ix_gac.run_LL(2*swap_days, LL_type='ix-gac')
        
        plt.figure(dpi=dpi)
        for comp in included_compounds:
            try: 
                p = plt.plot(results2['effluent results'].index,
                        results2['effluent results'][comp].values,
                        label=comp, #f"GAC-Lag: {comp}",
                        )
                
                plt.plot(results2['mid results'].index,
                        results2['mid results'][comp].values,
                        # label=f"IX-Lead: {comp}",
                        ls=':',
                        alpha=0.5,
                        color=p[0].get_color())
   
            except:
                pass
        
        plt.legend(loc='upper left', ncols=2, title="Solid: GAC-Lag, Dotted: IX-Lead")

        plt.title(f'IX Lead-GAC Lag (Swap Day: {2*swap_days})', fontsize=20)
        plt.ylabel('Concentration (ng/L)', fontsize=20)
        plt.xlabel('Time (days)', fontsize=20)
        plt.tight_layout()    
        
    
    
    
    