# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 09:02:55 2020

Example file outlining different options for PSDM - GAC modeling code

@author: UCChEJBB
"""

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

# =============================================================================
# Define the File to analyze
# can be in different folder, must be explicit for file location or can 
# use referential convention '../data/filename.xlsx' or 'C:/folders.../filename'
# =============================================================================
fn = 'Example_TCE.xlsx'

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

# =============================================================================
# Example run #1 - NO FOULING
# 
# chem_type = 'halogenated alkenes' #omitted, see Example 2 for how to use
# water_type = 'Organic Free'       #omitted, see Example 2 for how to use
# nr = radial collocation points
# nz = axial collocation points
# k_data = pandas dataframe containing K & 1/n data (ug/g)(L/ug)^(1/n) expected
# =============================================================================
print('Output files are saved to folder with example files')
print('Running Example 1\n','-'*50)
# File to write out results    
xlsx_fn1 = 'Example_TCE_' + carbons[0] + '_example1.xlsx'

for comp in compounds:
    print(comp)
    
    # SETTING UP PSDM simulation column
    column = PSDM.PSDM(column_info[carbons[0]], 
                   chem_data, 
                   raw_data,
                   nr=8,
                   nz=12, 
                   k_data=k_data,\
                   optimize=False
                   )
    
    print('EBCT: ', np.round(column.ebct, 2), ' min') 
    column.test_range = np.array([k_data[comp]['K']])
    column.xn_range = np.array([k_data[comp]['1/n']])
    
    # RUNNING the simulation
    #only results is used in this example
    _compound, _k, _xn, _ssqs, results = column.run_psdm_kfit(comp)
    
    ### plotting results 
    #units are converted from ug/L to mg/L 
    plt.plot(results.index, 
             results.values/1000., 
             label='PSDM')
    (raw_data[column.influent][comp]/1000.).plot.line(label='influent',linestyle=':')
    plt.legend()
    plt.title(comp+'\nK: '+repr(round(_k,2))+' ($\mu$g/g)(L/$\mu$g)$^{1/n}$ - 1/n: '+repr(round(_xn,3)))
    plt.xlabel(column_info[carbons[0]]['time'])
    plt.ylabel('Concentration - mg/L (ppm)')
    plt.savefig(comp+'_'+carbons[0]+'_example1.png', dpi=300)
    plt.close()
    
    tab2 = 'model-'+comp
    with pd.ExcelWriter(xlsx_fn1, engine='openpyxl') as writer:
        results.to_excel(writer, sheet_name=tab2)


# =============================================================================
# Example 1(a) - changing concentration of influent
# multiple EBCT
# =============================================================================
print('\nRunning Example 1a - Changing Concentration\nMultiple EBCT\n','-'*50)
raw_data2 = raw_data * 1.
raw_data2[column.influent] = 5000 # 5ppm, rather than 50 ppm in original file

xlsx_fn1a = 'Example_TCE_' + carbons[0] + '_example1a.xlsx'


for comp in compounds:
    print(comp)
    
    # SETTING UP PSDM simulation column
    column = PSDM.PSDM(column_info[carbons[0]], 
                   chem_data, 
                   raw_data2,
                   nr=8,
                   nz=12, 
                   k_data=k_data,\
                   optimize=False
                   )
    
    print('EBCT: ', np.round(column.ebct, 2), ' min') 
    column.test_range = np.array([k_data[comp]['K']])
    column.xn_range = np.array([k_data[comp]['1/n']])
    
    # RUNNING the simulation
    #only results is used in this example
    _compound, _k, _xn, _ssqs, results = column.run_psdm_kfit(comp)
    
    #Change Flowrate
    column_info2 = column_info[carbons[0]]  # makes copy to change flowrate
    column_info2['flrt'] = column_info2['flrt'] * 1.5 # internal units are mL/min, this just multiplies by 1.5
    column2 = PSDM.PSDM(column_info2, 
                   chem_data, 
                   raw_data2,
                   nr=8,
                   nz=12, 
                   k_data=k_data,\
                   optimize=False
                   )
    
    print('EBCT: ', np.round(column2.ebct, 2), ' min') 
    column2.test_range = np.array([k_data[comp]['K']])
    column2.xn_range = np.array([k_data[comp]['1/n']])
    _compound, _k, _xn, _ssqs, results2 = column2.run_psdm_kfit(comp)

    
    ### plotting results 
    plt.plot(results.index, 
             results.values/1000., 
             label='EBCT: '+repr(np.round(column.ebct,2)))
    plt.plot(results2.index, 
             results2.values/1000., 
             label='EBCT: '+repr(np.round(column2.ebct,2)))
    (raw_data2[column.influent][comp]/1000.).plot.line(label='influent',linestyle=':')
    plt.legend()
    plt.title(comp+'\nK: '+repr(round(_k,2))+' ($\mu$g/g)(L/$\mu$g)$^{1/n}$ - 1/n: '+repr(round(_xn,3)))
    plt.xlabel(column_info[carbons[0]]['time'])
    plt.ylabel('Concentration - mg/L (ppm)')
    plt.savefig(comp+'_'+carbons[0]+'_example1a.png', dpi=300)
    plt.close()
    
    tab2 = 'model-'+comp
    
    with pd.ExcelWriter(xlsx_fn1a, engine='openpyxl') as writer:
        results.to_excel(writer, sheet_name=tab2)




# =============================================================================
# Example 2 - variable influent
# variable influent data stored in "data_variable" sheet, must reload info.
# k data didn't change, so was reused from earlier
# =============================================================================
print('\nRunning Example 2 - Variable Influent\nWith Fouling\n','-'*50)

raw_data, column_info,\
compounds, carbons, = PSDM.process_input_file(fn,\
                                               data_sheet='data_variable',\
                                               column_sheet='columnSpecs'
                                              )

xlsx_fn2 = 'Example_TCE_' + carbons[0] + '_example2.xlsx'

chem_type = 'halogenated alkanes'
water_type = 'Rhine'
#note additional inputs in definition of column, below.

for comp in compounds:
    print(comp)
    
    # SETTING UP PSDM simulation column
    column = PSDM.PSDM(column_info[carbons[0]], 
                   chem_data, 
                   raw_data,
                   nr=8,
                   nz=12, 
                   k_data=k_data,
                   optimize=False
                   )
    
    print('EBCT: ', round(column.ebct, 2), ' min') 
    column.test_range = np.array([k_data[comp]['K']])
    column.xn_range = np.array([k_data[comp]['1/n']])
    
    # RUNNING the simulation
    #only results is used in this example
    _compound, _k, _xn, _ssqs, results = column.run_psdm_kfit(comp)
    
    column2 = PSDM.PSDM(column_info[carbons[0]], 
                   chem_data, 
                   raw_data,
                   nr=8,
                   nz=12, 
                   k_data=k_data,
                   chem_type=chem_type,
                   water_type=water_type,
                   optimize=False
                   )
    
    column2.test_range = np.array([k_data[comp]['K']])
    column2.xn_range = np.array([k_data[comp]['1/n']])
    
    # RUNNING the simulation
    #only results is used in this example
    _compound, _k, _xn, _ssqs, results2 = column2.run_psdm_kfit(comp)
    
    ### plotting results 
    #units are converted from ug/L to mg/L 
    plt.plot(results.index, 
             results.values/1000., 
             label='PSDM-no fouling')
    plt.plot(results2.index, 
             results2.values/1000., 
             label='PSDM-with fouling')
    (raw_data[column.influent][comp]/1000.).plot.line(label='influent',linestyle=':')
    plt.legend()
    plt.title(comp+'\nK: '+repr(round(_k,2))+' ($\mu$g/g)(L/$\mu$g)$^{1/n}$ - 1/n: '+repr(round(_xn,3)))
    plt.xlabel(column_info[carbons[0]]['time'])
    plt.ylabel('Concentration - mg/L (ppm)')
    plt.savefig(comp+'_'+carbons[0]+'_example2.png', dpi=300)
    plt.close()
    
    tab2 = 'model-'+comp
    with pd.ExcelWriter(xlsx_fn2, engine='openpyxl') as writer:

        results.to_excel(writer, sheet_name=tab2)



# =============================================================================
# Advanced Example 
# Using fitting capability
# =============================================================================
print('\nRunning Optimization Example (#3)\n', '-'*50, '\nWarning: May take a few minutes')

#reset files, just in case
chem_data = PSDM.process_input_data(fn, sheet_name='Properties') 
k_data = pd.read_excel(fn, sheet_name='Kdata',index_col=0) # K & 1/n data

raw_data, column_info,\
compounds, carbons, = PSDM.process_input_file(fn,\
                                               data_sheet='data_optimize',\
                                               column_sheet='columnSpecs'
                                              )
    
chem_type = 'halogenated alkanes'
water_type = 'Rhine'

for comp in compounds:
    print(comp)
    column = PSDM.PSDM(column_info[carbons[0]], 
                   chem_data, 
                   raw_data,
                   nr=8,
                   nz=12, 
                   # k_data=k_data,
                   chem_type=chem_type,
                   water_type=water_type,
                   optimize=False
                   )
    
    
    column.test_range = np.array([k_data[comp]['K']])
    column.xn_range = np.array([k_data[comp]['1/n']])
    
    # RUNNING the simulation
    #only results is used in this example
    _compound, _k, _xn, _ssqs, results = column.run_psdm_kfit(comp)
    
    # print(results['data'][np.linspace(0,175,36)])
    
    plt.plot(results.index, 
             results.values/1000,':',
             label='True Case: K: %i, 1/n: %.2f' % (k_data[comp]['K'],k_data[comp]['1/n'])
             )
    (raw_data[column_info[carbons[0]]['influentID']][comp]/1000).plot.line(label='influent',linestyle=':')
    (raw_data[carbons[0]][comp]/1000.).plot.line(label='effluent', marker='o',linestyle='--')
    
    
    column.xn_range = np.linspace(.3,.75,19) # increment 0.025
    # checks +/- 30% of estimated K (pm=30), 
    # and scans 1/n space in 0.01 step sizes (des_xn=0.01)
    column.run_all_smart(pm=30, des_xn=0.01)
    
    #Best fit was saved back to column.k_data
    print(column.k_data)
    
    column.xn_range = np.array([column.k_data[comp]['1/n']])
    column.test_range = np.array([column.k_data[comp]['K']])
    _compound, _k, _xn, _ssqs, results2 = column.run_psdm_kfit(comp)
    
    plt.plot(results2.index,
             results2.values/1000.,
             label='Best Fit')
    plt.legend()
    plt.title(comp+'\nK: '+repr(round(_k,2))+' ($\mu$g/g)(L/$\mu$g)$^{1/n}$ - 1/n: '+repr(round(_xn,3)))
    plt.xlabel(column_info[carbons[0]]['time'])
    plt.ylabel('Concentration - mg/L (ppm)')
    plt.savefig(comp+'_'+carbons[0]+'_example3.png', dpi=300)
    plt.close()
    
    '''
# =============================================================================
# ALTERNATE OPTIMIZER
# Depending on selections this can take 1+ hours per compound !!!!
# Replaces lines 306-310 above
# Commented out so it does not run initially
# =============================================================================
    '''
    # column.xn_range = np.linspace(0.3, 0.75, 46) 
    # '''46 = interval 0.01 (19 is 0.025, for 0.3-0.75)'''
    # column.test_range = np.linspace(1, 5, 41)
    # ''' 41 gives you 0.1 increment on K_multiplier
    # K would not be expected to be lower than 1, but upper bound of 3 to 4 may 
    # be sufficient to capture any real values'''
    # column.run_all()
    # column.run_all(optimize='brute') #run only this line or the line above, not both
    # print(column.k_data) # best fits are stored back to column.k_data
    # '''
    # run_all() by default uses a multi-staged brute force approach
    # This subdivides the space into coarse grids, and performs a series of 
    # refining steps. 'optimizer=brute' searches every combination of 
    # xn_range and test_range, while the default behavior omits sections that
    # are not near the minimum. 
    
    # 'optimize=brute' could be used after run_all_smart() or the default to 
    # ensure that the K & 1/n values around a minimum were thoroughly searched.
    
    # Performance: brute (1+ hours per compound)
    #              staged (~5-10+ minutes per compound)
    #              using run_all_smart (~2-15 minutes per compound)
    # '''

# =============================================================================
# Final Notes: If you got this far, you may be interested in using this for 
# larger datasets. In the above examples, each of the examples is set up
# for looped execution. If you provide datasets with more than one compound
# in the data sheet (and associated properties) then this will loop over each
# compound individually
# =============================================================================
