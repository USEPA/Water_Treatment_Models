# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:40:42 2019

@author: JBurkhar
"""

#test_PSDM.py

import numpy as np
import pandas as pd
import PSDM
import matplotlib.pyplot as plt
plt.close()
from scipy import interpolate

# read in test input and output data 
# output tests were run using AdDesignS, with the included .dat file
input_data = pd.read_excel('test.xlsx',sheet_name = 'data',index_col = [0])

output5 = pd.read_excel('test.xlsx', sheet_name = '5', index_col = [0], header = [0])
Cout_5 = output5['C, mg/L']
Cout_5 = Cout_5[Cout_5>0.]
plt.plot(Cout_5.index.values/1440, Cout_5.values)
Cout_5f = interpolate.interp1d(Cout_5.index.values/1440, Cout_5.values)

output25 = pd.read_excel('test.xlsx', sheet_name = '25', index_col = [0], header = [0])
Cout_25 = output25['C, mg/L']
Cout_25 = Cout_25[Cout_25>0.]
plt.plot(Cout_25.index.values/1440, Cout_25.values)
Cout_25f = interpolate.interp1d(Cout_25.index.values/1440, Cout_25.values)

output50 = pd.read_excel('test.xlsx', sheet_name = '50', index_col = [0], header = [0])
Cout_50 = output50['C, mg/L']
Cout_50 = Cout_50[Cout_50>0.] 
plt.plot(Cout_50.index.values/1440, Cout_50.values)
Cout_50f = interpolate.interp1d(Cout_50.index.values/1440, Cout_50.values)

plot_line = np.arange(np.max(Cout_5.index.values/1440))

#transform input data into appropriate data structures
column_prop = input_data.transpose()[['epor','diam', 'L','wt','rhop','rhof','flrt','rad','influentID']].transpose()
column_prop.loc['tortu'] = 1
column_prop.loc['psdfr'] = 5
test_data = input_data.transpose()[['MW','MolarVol']].transpose()
data_store = input_data.transpose()[['K','kf','dp','1/n']].transpose()
data_store.loc['q'] = 1
data_store.loc['brk'] = 100
data_store.columns = pd.MultiIndex.from_tuples([('Test','Test')])

xn = data_store['Test']['Test']['1/n']




#create influent data structures for 5 & 50 mg/L in ug/L (x1000)
idx = pd.MultiIndex.from_tuples([('concentration','Test','Test')]) #matches function's expectations
time = np.arange(175)
test5 = pd.DataFrame(5.*1e6, index = time, columns = idx) # 5ppm in ppt
test25 = pd.DataFrame(25.*1e6, index = time, columns = idx) # 25ppm in ppt
test50 = pd.DataFrame(50.*1e6, index = time, columns = idx) # 50ppm in ppt


test_column = PSDM.PSDM(column_prop['Test'],\
                        test_data,\
                        test5['concentration'],\
                        xn=xn,\
                        conc_type='ng',\
                        chem_type='chlorinated hydrocarbon',\
                        water_type='Rhine',\
                        k_data=data_store['Test'],\
                        test_range=data_store.loc['K'].values,\
                        xn_range=data_store.loc['1/n'].values,\
                        optimize=False)

comp, K, xn, ssq, md = test_column.run_psdm_kfit('Test')
plt.plot(md.index.values, md.values/1e6,'--')
md_5f = interpolate.interp1d(md.index.values, md.values.flatten())

idxs = Cout_5.index.values/1440
ssq5 = ((Cout_5f(idxs) - md_5f(idxs)/1e6)**2).sum()
print('Rhine River Test')
print('5: ' + repr(ssq5))


test_column = PSDM.PSDM(column_prop['Test'],\
                        test_data,\
                        test25['concentration'],\
                        xn=xn,\
                        conc_type='ng',\
                        chem_type='chlorinated hydrocarbon',\
                        water_type='Rhine',\
                        k_data=data_store['Test'],\
                        test_range=data_store.loc['K'].values,\
                        xn_range=data_store.loc['1/n'].values,\
                        optimize=False)

comp, K, xn, ssq, md = test_column.run_psdm_kfit('Test')
plt.plot(md.index.values, md.values/1e6,'--')
md_25f = interpolate.interp1d(md.index.values, md.values.flatten())

idxs = Cout_25.index.values/1440
ssq25 = ((Cout_25f(idxs) - md_25f(idxs)/1e6)**2).sum()
print('25: ' + repr(ssq25))


test_column = PSDM.PSDM(column_prop['Test'],\
                        test_data,\
                        test50['concentration'],\
                        xn=xn,\
                        conc_type='ng',\
                        chem_type='chlorinated hydrocarbon',\
                        water_type='Rhine',\
                        k_data=data_store['Test'],\
                        test_range=data_store.loc['K'].values,\
                        xn_range=data_store.loc['1/n'].values,\
                        optimize=False)

comp, K, xn, ssq, md = test_column.run_psdm_kfit('Test')
plt.plot(md.index.values, md.values/1e6,'--')
md_50f = interpolate.interp1d(md.index.values, md.values.flatten())

idxs = Cout_50.index.values/1440
ssq50 = ((Cout_50f(idxs) - md_50f(idxs)/1e6)**2).sum()
print('50: ' + repr(ssq50))



print('No Fouling Test')
       
# =============================================================================
# no fouling considered
# =============================================================================
    
output5 = pd.read_excel('test.xlsx', sheet_name = '5no', index_col = [0], header = [0])
Cout_5 = output5['C, mg/L']
Cout_5 = Cout_5[Cout_5>0.]
plt.plot(Cout_5.index.values/1440, Cout_5.values)
Cout_5f = interpolate.interp1d(Cout_5.index.values/1440, Cout_5.values)

output25 = pd.read_excel('test.xlsx', sheet_name = '25no', index_col = [0], header = [0])
Cout_25 = output25['C, mg/L']
Cout_25 = Cout_25[Cout_25>0.]
plt.plot(Cout_25.index.values/1440, Cout_25.values)
Cout_25f = interpolate.interp1d(Cout_25.index.values/1440, Cout_25.values)

output50 = pd.read_excel('test.xlsx', sheet_name = '50no', index_col = [0], header = [0])
Cout_50 = output50['C, mg/L']
Cout_50 = Cout_50[Cout_50>0.] 
plt.plot(Cout_50.index.values/1440, Cout_50.values)
Cout_50f = interpolate.interp1d(Cout_50.index.values/1440, Cout_50.values)

test_column = PSDM.PSDM(column_prop['Test'],\
                        test_data,\
                        test5['concentration'],\
                        xn=xn,\
                        conc_type='ng',\
                        k_data=data_store['Test'],\
                        test_range=data_store.loc['K'].values,\
                        xn_range=data_store.loc['1/n'].values,\
                        optimize=False)

comp, K, xn, ssq, md = test_column.run_psdm_kfit('Test')
plt.plot(md.index.values, md.values/1e6,'--')
md_5f = interpolate.interp1d(md.index.values, md.values.flatten())

idxs = Cout_5.index.values/1440
ssq5 = ((Cout_5f(idxs) - md_5f(idxs)/1e6)**2).sum()
print('5: ' + repr(ssq5))


test_column = PSDM.PSDM(column_prop['Test'],\
                        test_data,\
                        test25['concentration'],\
                        xn=xn,\
                        conc_type='ng',\
                        k_data=data_store['Test'],\
                        test_range=data_store.loc['K'].values,\
                        xn_range=data_store.loc['1/n'].values,\
                        optimize=False)

comp, K, xn, ssq, md = test_column.run_psdm_kfit('Test')
plt.plot(md.index.values, md.values/1e6,'--')
md_25f = interpolate.interp1d(md.index.values, md.values.flatten())

idxs = Cout_25.index.values/1440
ssq25 = ((Cout_25f(idxs) - md_25f(idxs)/1e6)**2).sum()
print('25: ' + repr(ssq25))


test_column = PSDM.PSDM(column_prop['Test'],\
                        test_data,\
                        test50['concentration'],\
                        xn=xn,\
                        conc_type='ng',\
                        k_data=data_store['Test'],\
                        test_range=data_store.loc['K'].values,\
                        xn_range=data_store.loc['1/n'].values,\
                        optimize=False)

comp, K, xn, ssq, md = test_column.run_psdm_kfit('Test')
plt.plot(md.index.values, md.values/1e6,'--')
md_50f = interpolate.interp1d(md.index.values, md.values.flatten())

idxs = Cout_50.index.values/1440
ssq50 = ((Cout_50f(idxs) - md_50f(idxs)/1e6)**2).sum()
print('50: ' + repr(ssq50))

print('Keweenaw Waterway (Portage) Test')
# =============================================================================
# Keweenaw Waterway (Portage) Test
# =============================================================================
output50 = pd.read_excel('test.xlsx', sheet_name = '50K', index_col = [0], header = [0])
Cout_50 = output50['C, mg/L']
Cout_50 = Cout_50[Cout_50>0.] 
plt.plot(Cout_50.index.values/1440, Cout_50.values)
Cout_50f = interpolate.interp1d(Cout_50.index.values/1440, Cout_50.values)

test_column = PSDM.PSDM(column_prop['Test'],\
                        test_data,\
                        test50['concentration'],\
                        xn=xn,\
                        conc_type='ng',\
                        chem_type='chlorinated hydrocarbon',\
                        water_type='Portage',\
                        k_data=data_store['Test'],\
                        test_range=data_store.loc['K'].values,\
                        xn_range=data_store.loc['1/n'].values,\
                        optimize=False)

comp, K, xn, ssq, md = test_column.run_psdm_kfit('Test')
plt.plot(md.index.values, md.values/1e6,'--')
md_50f = interpolate.interp1d(md.index.values, md.values.flatten())

idxs = Cout_50.index.values/1440
ssq50 = ((Cout_50f(idxs) - md_50f(idxs)/1e6)**2).sum()
print('50: ' + repr(ssq50))

plt.xlim((0, max(Cout_5.index.values/1440)))
plt.ylim((0,70))