# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 13:18:37 2019

@author: UCChEJBB
"""


import PSDM_tools
import pandas as pd

#Example data, user can save this data to Excel file for structure, or provide 
# data from Excel file for real use
# 
# Currently no unit-handling is performed, units output will be related to input

isotherm_data = pd.DataFrame([[3.37, 3.37, 0.000,'carbon1','test1'],
                              [3.37, 3.27, 0.001,'carbon1','test1'],
                              [3.37, 2.77, 0.010,'carbon1','test1'],
                              [3.37, 1.86, 0.100,'carbon1','test1'],
                              [3.37, 1.33, 0.500,'carbon1','test1'],
                              [100, 100, 0.000,'carbon2','test2'],
                              [100, 33, 168,'carbon2','test2'],
                              [100, 20.5, 265,'carbon2','test2'],
                              [100, 12.5, 398,'carbon2','test2'],
                              [100, 8.0, 548,'carbon2','test2'],
                              [100, 5.2, 740,'carbon2','test2'],
                              [100, 3.5, 965,'carbon2','test2'],
                              [23.08, 23.08, 0.000,'carbon1','test2'],
                              [23.08, 20.3, 28.71,'carbon1','test2'],
                              [23.08, 15.9, 84.53,'carbon1','test2'],
                              [23.08, 11.7, 183.63,'carbon1','test2'],
                              [23.08, 7.4, 340.97,'carbon1','test2'],
                              [23.08, 4.2, 651.2,'carbon1','test2'],
                              [23.08, 1.8, 1934.98,'carbon1','test2']],\
                              columns = ['C0','Ce','mass','carbon','compound'])  

compounds = set(isotherm_data['compound'])
carbons = set(isotherm_data['carbon'])
# print(compounds, carbons)

# Example with input that contains both possible carbons and "compounds", lited as 
for compound in compounds:
    for carbon in carbons:
        sub_data = isotherm_data[((isotherm_data['carbon']==carbon) &\
                                  (isotherm_data['compound']==compound))]
        if len(sub_data.index) > 0:
            # PSDM_tools.isotherm_fit(sub_data, isotherm='freundlich')
            # PSDM_tools.isotherm_fit(sub_data, isotherm='langmuir')
            # PSDM_tools.isotherm_fit(sub_data, isotherm='RedlichPeterson')
            pass
        
# =============================================================================
# Example Chen, Information 2015
# =============================================================================

isotherm_data2 = pd.DataFrame([[0.49, 0.021, 'm1'],
                               [0.78, 0.39, 'm1'],
                               [3.67, 1.91, 'm1'],
                               [9.03, 3.88, 'm1'],
                               [23.98, 4.82, 'm1'],
                               [53.25, 5.01, 'm1'], 
                               [0.061, 0.93, 'm2'],
                               [0.023, 1.91, 'm2'],
                               [0.35, 8.56, 'm2'],
                               [2.45, 17.18, 'm2'],
                               [16.43, 19.74, 'm2'],
                               [44.40, 20.91, 'm2'],
                               [0.037, 0.98, 'm3'],
                               [0.013, 1.92, 'm3'],
                               [0.036, 9.26, 'm3'],
                               [0.15, 19.49, 'm3'],
                               [5.74, 42.22, 'm3'],
                               [29.5, 44.3, 'm3'], 
                               [0.052, 0.95, 'm4'],
                               [0.017, 1.94, 'm4'],
                               [0.034, 9.26, 'm4'],
                               [0.33, 18.97, 'm4'],
                               [4.12, 44.95, 'm4'],
                               [26.09, 52.14, 'm4']],\
                              columns=['Ce','q','case'])

cases = set(isotherm_data2['case'])    
for case in cases:
    print(case)
    sub_data = isotherm_data2[isotherm_data2['case']==case]
    PSDM_tools.isotherm_fit(sub_data, isotherm='freundlich')
    PSDM_tools.isotherm_fit(sub_data, isotherm='langmuir')

    