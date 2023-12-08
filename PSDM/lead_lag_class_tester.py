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






##
print('\nTesting\n')
print('x'*30,'\n')


test_gac_gac = False
test_ix_ix = False
test_gac_ix = False
test_ix_gac = True

### Generic Testing
####   TEST GAC
if test_gac_gac:
    
    
    fn = srt_dir + '/examples/Example_TCE.xlsx'

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

    raw_data /= 10 # reduce overall concentration
        
        
    swap_days = 133
     
    plt.figure()
    
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
                 results.values/1e3,
                 ls=':',
                 label=f"{comp} - Single")
        
        
    # =============================================================================
    # run lead-lag for gac    
    # =============================================================================
    LL_gac = LLobj(raw_data=raw_data,
                   properties=chem_data,
                   column_info=column_info[carbons[0]],
                   k_data=k_data,
                   )
    
    results_dict = LL_gac.run_LL(swap_days, LL_type='gac')
    for comp in compounds:
        
        mid = results_dict['mid results'][comp]
        plt.plot(mid.index, mid.values/1e3, label=f"{comp} - Mid")
        
        effl = results_dict['effluent results'][comp]
        plt.plot(effl.index, effl.values/1e3, label=f"{comp} - Effl")
        
    
    plt.legend()
    plt.title(f"Swap Days: {swap_days}")
    plt.ylabel('Concentration (mg/L)')
    plt.xlabel('Time (days)')
    plt.xlim((0, swap_days))

### Test IX
if test_ix_ix: 
    ix_fn = 'test_iex/lag_iex.xlsx'
    
    LL_ix = LLobj(ix_filename=ix_fn, ix_nr=10, ix_nz=15)
    
    swap_days = 500
    results = LL_ix.run_LL(swap_days, LL_type='ix')
    print(results)
    
    for comp in results['effluent'].keys():
        temp_res = results['effluent results'][comp]
        if comp not in ['CHLORIDE','SULFATE','BICARBONATE','NITRATE']:
            plt.plot(temp_res.index,
                     temp_res.values, ## file had mg/L specified for PFAS
                     label=comp)
    plt.legend()
    plt.title(f"Swap Days: {swap_days}")
    plt.ylabel('Concentration (ng/L)')
    plt.xlabel('Time (days)')
    
    plt.figure()
    
    for comp in results['effluent'].keys():
        temp_res = results['effluent results'][comp]
        if comp in ['CHLORIDE','SULFATE','BICARBONATE','NITRATE']:
            mw = LL_ix.ix_data['ions'][LL_ix.ix_data['ions']['name'] == comp]['mw'].values[0]
            val = LL_ix.ix_data['ions'][LL_ix.ix_data['ions']['name'] == comp]['valence'].values[0]
            
            plt.plot(temp_res.index,
                     temp_res.values * mw / val,
                     label=comp)
    
    plt.legend()
    plt.xlim((0,4))
    plt.title(f"Swap Days: {swap_days}")
    plt.ylabel('Concentration (mg/L)')
    plt.xlabel('Time (days)')
    
    
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
    if test_gac_ix:
        LL_gac_ix = LLobj(raw_data=raw_data,
                            properties=chem_data,
                            column_info=column_info[carbons[0]],
                            k_data=k_data,
                            ix_filename=ix_fn
                        )
        
        results = LL_gac_ix.run_LL(500, LL_type='gac-ix')
    
        print(results)
        results['effluent results'].plot.line()
    
    
    
    ### IX-GAC 
    if test_ix_gac:
        LL_ix_gac = LLobj(raw_data=raw_data,
                            properties=chem_data,
                            column_info=column_info[carbons[0]],
                            k_data=k_data,
                            ix_filename=ix_fn,
                            water_type='Rhine',
                            chem_type='PFAS',
                        )
        results2 = LL_ix_gac.run_LL(1500, LL_type='ix-gac')
        
        
        plt.figure()
        for comp in results2['effluent results'].columns:#['PFHxA']:#'PFHpA', 'PFHxS', 'PFBA']:
            # plt.figure()
            plt.plot(results2['mid results'].index,
                     results2['mid results'][comp].values,
                     label=f"IX-Lead: {comp}",
                     ls=':')
            plt.plot(results2['effluent results'].index,
                     results2['effluent results'][comp].values,
                     label=f"GAC-Lag: {comp}")
        
        plt.legend()

        plt.title('IX Lead-GAC Lag')
        plt.ylabel('Concentration (ng/L)')
        plt.xlabel('Time (days)')
        plt.tight_layout()    
        
    
    
    
    