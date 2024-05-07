# -*- coding: utf-8 -*-
"""

@author: UCChEJBB - Jonathan Burkhardt

Contains convenience functions and variables used by PSDM. 
Moved in order to keep PSDM.py file size shorter, and better manage code.

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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
from scipy import special
from scipy.interpolate import interp1d
from scipy import sparse

lpg = 3.785411784 # liter per gallon conversion
cm_per_ft = 2.54 * 12
min_per_day = 24 * 60
gm_per_lb = 453.59237

# =============================================================================
# Input reader functions
# =============================================================================
def process_input_file(filename, data_sheet='data',\
                          column_sheet='column info', rept_lim=0):
    '''
    Input file reader. Processes data into formats/structures used by PSDM.

    Parameters
    ----------
    filename : string
        Example: 'file.xlsx'
        filename of file containing input parameters and information. 
    data_sheet : string, optional
        Sheet name where influent and effluent data is found. 
        The default is 'data'.
    column_sheet : string, optional
        Sheet name where the column data can be found. 
        The default is 'column info'.
    rept_lim : float/int, optional
        Can be used to filter out values bellow a reporting limit. 
        The default is 0.

    Returns
    -------
    rawdata_df : Pandas Dataframe
        DESCRIPTION.
    column_data : TYPE
        DESCRIPTION.
    compounds : TYPE
        DESCRIPTION.
    carbons : TYPE
        DESCRIPTION.

    '''
    
    data_df = pd.read_excel(filename, sheet_name=data_sheet, index_col=0)
    c_data = pd.read_excel(filename, sheet_name=column_sheet, index_col=0)
    
    
    column_dict = {'rad': 0, 'flrt': 0, 'epor': 0,\
                   'psdfr': 5, 'rhop': 0, 'rhof': 0, 'L': 0, \
                   'wt': 0, 'diam': 0, 'tortu':1, 'influentID':'influent',\
                   'effluentID': 'effluent', 'units':'ng', 'time':'days',\
                   'mass_mul': 1., 't_mult': 1., 'flow_mult': 1., \
                   'flow_type': 'ml'}
    column_data = pd.DataFrame.from_dict(column_dict, orient='index')
    column_data.columns = [c_data['value']['carbonID']]
    carbons=list(column_data.columns)

    # recalculate required values to default units
    #dimensions in cm (pardicle radius (rad), length (L) and diameter (diam) )
    dims = {'L': 'length', 'diam': 'diameter', 'rad': 'radius'}
    for dim in dims.keys():
        val = dims[dim]
        unit = (c_data['units'][val]).lower() # converts everything to lowercase
        tmp_value = c_data['value'][val]
        if unit == 'cm':
            column_data.loc[dim] = tmp_value * 1.
        elif unit == 'mm':
            column_data.loc[dim] = tmp_value / 10.
        elif unit == 'in':
            column_data.loc[dim] = tmp_value * 2.54
        elif unit == 'm':
            column_data.loc[dim] = tmp_value * 100.
        elif unit == 'ft':
            column_data.loc[dim] = tmp_value * cm_per_ft
    # densities (rhop and rhof) = g/ml
    rhos = {'rhop': 'particleDensity', 'rhof': 'apparentDensity'}
    for rho in rhos.keys():
        val = rhos[rho]
        unit = (c_data['units'][val]).lower()
        tmp_value = c_data['value'][val] 
        if unit == 'g/ml' or unit == 'g/cm3':
            column_data.loc[rho] = tmp_value * 1.
        elif unit == 'lb/ft3':
            column_data.loc[rho] = tmp_value / 62.427973725314 #magic #
    #flowrate = ml/min
    unit = (c_data['units']['flowrate']).lower()
    flrt = c_data['value']['flowrate']
    if unit == 'ml/min' or unit == 'cm3/min':
        column_data.loc['flrt'] = flrt * 1.
    elif unit == 'lpm' or unit == 'l/min':
        column_data.loc['flrt'] = flrt * 1e3
    elif unit == 'gpm' or unit == 'gal/min':
        column_data.loc['flrt'] = flrt * lpg * 1e3
    elif unit == 'mgd':
        column_data.loc['flrt'] = flrt * lpg * 1e3 * 1e6 / min_per_day # 1e6 = million, 1e3 L-> ml
    elif unit == 'm3/min':
        column_data.loc['flrt'] = flrt * 1e6
    elif unit == 'ft3/min':
        column_data.loc['flrt'] = flrt * cm_per_ft**3
    # carbon mass = grams
    unit = (c_data['units']['weight']).lower()
    wt = c_data['value']['weight']
    if unit == 'g' or unit == 'gm':
        column_data.loc['wt'] = wt * 1.
    elif unit == 'lb':
        column_data.loc['wt'] = wt * gm_per_lb
    elif unit == 'kg':
        column_data.loc['wt'] = wt * 1e3
    elif unit == 'ton':
        column_data.loc['wt'] = wt * 2e3 * gm_per_lb #2000 lbs per ton
    elif unit == 'mton' or unit == 'metricton':
        column_data.loc['wt'] = wt * 1e3 * 1e3 #1000 kg per metric ton
    #unitless
    column_data.loc['tortu'] = c_data['value']['tortuosity'] * 1.
    column_data.loc['epor'] = c_data['value']['porosity'] * 1.
    column_data.loc['psdfr'] = c_data['value']['psdfr'] * 1.
    column_data.loc['influentID'] = c_data['value']['influentID']
    column_data.loc['effluentID'] = c_data['value']['effluentID']
    #add more units handling # some redundancy from __init__ in psdm
    units = (c_data['value']['units']).lower() #may need to add error handling for if it is under units
    if units == 'ug' or units == 'ng':
        column_data.loc['units'] = units
    elif units == 'ppb' or units == 'ug/l':
        column_data.loc['units'] = 'ug'
    elif units == 'ppt' or units == 'ng/l':
        column_data.loc['units'] = 'ng'
    
    if column_data[carbons[0]]['units'] == 'ug':
        column_data.loc['mass_mul'] = 1.
    elif column_data[carbons[0]]['units'] == 'ng':
        column_data.loc['mass_mul'] = 1e-3    
            
    #process data_df into rawdata_df
    compounds = list(set(data_df['compound'].values))
    types = list(set(data_df.index.values))
    influent = column_data[column_data.columns[0]]['influentID']
    effluent = column_data[column_data.columns[0]]['effluentID']
    if influent not in types:
        print('Error: No data associated with influent found')
    if effluent not in types:
        print('Error: No data associated with effluent found')
    
    #moves some of the content from __init__ of PSDM
    if 'time' in c_data.index:
        if c_data['value']['time'] == 'days' or \
           c_data['value']['time'] == 'day':
            column_data.loc['time'] = 'days'
            column_data.loc['t_mult'] = min_per_day
        elif c_data['value']['time'] == 'hour' or \
             c_data['value']['time'] == 'hr' or \
             c_data['value']['time'] == 'hrs':
            column_data.loc['time'] = 'hour'
            column_data.loc['t_mult'] = 60. #min_per_hour
        elif c_data['value']['time'] == 'min':
            column_data.loc['time'] = 'min'
            column_data.loc['t_mult'] = 1. # minutes
    else: # assume days
        column_data.loc['time'] = 'days'
        column_data.loc['t_mult'] = min_per_day
            
    #previous code converts flowrate into ml
    column_data.loc['flow_mult'] = 1e-3
        
    #Restructure Data and Clean
    data_df.columns = [i.lower() for i in data_df.columns]
    if 'time' in data_df.columns:
        a = data_df.pivot_table(index=['time'],columns=['type','compound'])
    elif 'days' in data_df.columns:
        a = data_df.pivot_table(index=['days'],columns=['type','compound'])
    a[a < rept_lim] = 0. # assume zero below RepLim, LoQ, LoR
    rawdata_df = a['concentration']
    
    infl_idx = list(rawdata_df[influent].dropna().index)
    effl_idx = list(rawdata_df[effluent].dropna().index)
    
    infl_dict = {}
    effl_dict = {}
    for comp in compounds:
        infl_dict[comp] = interp1d(infl_idx,\
                                    rawdata_df[influent][comp][infl_idx].values, \
                                    fill_value = 'extrapolate')
        effl_dict[comp] = interp1d(effl_idx,\
                                    rawdata_df[effluent][comp][effl_idx].values, \
                                    fill_value = 'extrapolate')
    #handles missing initial point
    #should be handled in Excel file, but just in case
    if not 0 in infl_idx:
        infl_idx = np.sort(infl_idx + [0])
    if not 0 in effl_idx:
        effl_idx = np.sort(effl_idx + [0])
    idx = list(rawdata_df.index)
    if not 0 in idx:
        idx = np.sort(idx + [0])
        rawdata_df.loc[0] = rawdata_df.loc[infl_idx[-1]] * 0.

    diff_infl = list(set(idx).difference(infl_idx))
    diff_effl = list(set(idx).difference(effl_idx))
    
    #don't need to worry about values if day index is less than 0
    diff_infl = [i for i in diff_infl if i >= 0 ]
    diff_effl = [i for i in diff_effl if i >= 0 ]
    
    #fill in missing data with interpolation from known data
    #could also just replace with average of known data?
    if len(diff_infl) > 0:
        for comp in compounds:
            for i in diff_infl:
                rawdata_df[influent][comp][i] = infl_dict[comp](i)
    if len(diff_effl) > 0:
        for comp in compounds:
            for i in diff_effl:
                rawdata_df[effluent][comp][i] = effl_dict[comp](i)
    
    #resort index to makes sure things are in order
    rawdata_df = rawdata_df.sort_index()
    rawdata_df[influent][rawdata_df[influent]==0.] = 0.001 #makes this value very small, but non-zero
    rawdata_df = rawdata_df[[influent,effluent]]
    
    top_lvl = rawdata_df.columns.levels[0]
    new_lvl = []
    for i in top_lvl:
        if i == effluent:
            new_lvl.append(carbons[0])
        else:
            new_lvl.append(i)

    rawdata_df.columns = rawdata_df.columns.set_levels(new_lvl, level=0)

    return rawdata_df, column_data, compounds, carbons

def process_column_data(filename):
    '''
    Parameters
    ----------
    filename : string, example = 'file.xlsx'
        Filename for Excel file containing column information. 
        Assumes Sheet1 contains relevant information.
        
        Maintained for older code that used this function:
            use 'process_input_file' instead. It integrates column and raw 
            data handling into one file

    Returns
    -------
    column_data :  pandas dataframe
        column dataframe used by PSDM code.
        
        if error, returns emmpty dataframe.

    '''
    try:
        column_data = pd.read_excel(filename, sheet_name='Sheet1', index_col=0)
    except:
        column_data = pd.DataFrame()
    return column_data

def process_input_data(filename, sheet_name='Sheet1'):
    '''
    Parameters
    ----------
    filename : string, example = 'file.xlsx'
        filename for Excel file that contains table of compound parameters.
        Molecular weight, molar volume, and ..... are used in PSDM functions.
    sheet_name : string, example: sheet_name='Sheet1', optional
        sheet name of sheet with required data, The default is 'Sheet1'      
        
    Returns
    -------
    comp_data : pandas dataframe
        compound information dataframe.
        
        if error, returns empty dataframe.

    '''
    try:
        comp_data = pd.read_excel(filename,sheet_name=sheet_name, index_col=0)
    except:
        comp_data = pd.DataFrame()
    return comp_data

def process_raw_data(filename, rept_lim=0.):
    '''
    USE 'process_input_file' instead

    Parameters
    ----------
    filename : string, 
        File containing raw data to be used in analysis.
    rept_lim : float, optional
        Reporting limit for analysis. This value can be used to filter 
        analytical results below a limit if desired. The default is 0.

    Returns
    -------
    rawdata_df : TYPE
        DESCRIPTION.
    compounds : TYPE
        DESCRIPTION.
    carbons : TYPE
        DESCRIPTION.

    '''
    data_df = pd.read_csv(filename) 
    # compounds and carbons to consider     
    compounds = list(set(data_df['compound'].values))
    carbons = list(set(data_df['type'].values))
    influents = [i for i in carbons if 'Influent' in i]
    carbons = [i for i in carbons if 'Influent' not in i] #strip out the 'Influent' names
    
    #Restructure Data and Clean
    a = data_df.pivot_table(index = ['time'],columns = ['type','compound'])
    a = a.fillna(0.)
    a[a<=rept_lim] = 0. # min reporting limit = 10.3, assume zero
    rawdata_df = a['concentration']
    rawdata_df[influents][rawdata_df[influents]==0.] = 0.001 #makes this value very small, but non-zero
    
    return rawdata_df, compounds, carbons

def process_raw_data_new(filename, rept_lim=0):
    '''
    USE 'process_input_file' instead

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    rept_lim : TYPE, optional
        DESCRIPTION. The default is 0.

    Returns
    -------
    rawdata_df : TYPE
        DESCRIPTION.
    column_data : TYPE
        DESCRIPTION.
    compounds : TYPE
        DESCRIPTION.
    carbons : TYPE
        DESCRIPTION.

    '''
    data_df = pd.read_excel(filename, sheet_name = 'data', index_col = 0)
    column_data = pd.read_excel(filename, sheet_name = 'column info', index_col = 0)
    
    compounds = list(set(data_df['compound'].values))
    carbons = column_data.columns
    influents = column_data[carbons[0]]['influentID']
    try:
        effluents = column_data[carbons[0]]['effluentID']
    except:
        effluents = 'Intermediate'
    
    #Restructure Data and Clean
    a = data_df.pivot_table(index = ['time'],columns = ['type','compound'])
    a[a<=rept_lim] = 0. # min reporting limit = 10.3, assume zero
    rawdata_df = a['concentration']
    
    infl_idx = list(rawdata_df[influents].dropna().index)
    effl_idx = list(rawdata_df[effluents].dropna().index)
    
    infl_dict = {}
    effl_dict = {}
    for comp in compounds:
        infl_dict[comp] = interp1d(infl_idx,\
                                    rawdata_df[influents][comp][infl_idx].values, \
                                    fill_value = 'extrapolate')
        effl_dict[comp] = interp1d(effl_idx,\
                                    rawdata_df[effluents][comp][effl_idx].values, \
                                    fill_value = 'extrapolate')
    
    if not 0 in infl_idx:
        infl_idx = np.sort(infl_idx + [0])
    if not 0 in effl_idx:
        effl_idx = np.sort(effl_idx + [0])
    idx = list(rawdata_df.index)
    if not 0 in idx:
        idx = np.sort(idx + [0])
        rawdata_df.loc[0] = rawdata_df.loc[infl_idx[-1]] * 0.
        
    diff_infl = list(set(idx).difference(infl_idx))
    diff_effl = list(set(idx).difference(effl_idx))
    
    #don't need to worry about values if day index is less than 0
    diff_infl = [i for i in diff_infl if i >= 0 ]
    diff_effl = [i for i in diff_effl if i >= 0 ]
    
    #fill in missing data with interpolation from known data
    #could also just replace with average of known data?
    if len(diff_infl) > 0:
        for comp in compounds:
            for i in diff_infl:
                rawdata_df[influents][comp][i] = infl_dict[comp](i)
    if len(diff_effl) > 0:
        for comp in compounds:
            for i in diff_effl:
                rawdata_df[effluents][comp][i] = effl_dict[comp](i)
    
    #resort index to makes sure things are in order
    rawdata_df = rawdata_df.sort_index()
    
    rawdata_df[influents][rawdata_df[influents]==0.] = 0.001 #makes this value very small, but non-zero
    
    rawdata_df = rawdata_df[[influents,effluents]]
    
    top_lvl = rawdata_df.columns.levels[0]
    new_lvl = []
    for i in top_lvl:
        if i == effluents:
            new_lvl.append(carbons[0])
        else:
            new_lvl.append(i)
    
    
    rawdata_df.columns = rawdata_df.columns.set_levels(new_lvl, level=0)
    
    return rawdata_df, column_data, compounds, carbons


# =============================================================================
# Fouling Parameters
# =============================================================================
# edited to match AdDesignS 1/min, not 1/day
# 
# PFAS Specific Data presented in Burkhardt et al. 2020, submitted to 
# AWWA Water Science - set to halogenated alkanes until published
foul_params = {'water':{'Organic Free':[1.,          0.,    0.,       0.],
                              'Rhine': [0.35,  -6.15e-8,  0.65, -8.93e-5],
                            'Portage': [0.510, -9.21e-7, 0.490, -2.80e-5], 
                          'Karlsruhe': [0.65,  -6.71e-7,  0.35, -1.00e-4], 
                             'Wausau': [0.83,  -9.12e-7,  0.17, -2.65e-4], 
                           'Houghton': [0.66, - 1.55e-7,  0.34, -7.29e-5]}, 
               'chemical':{'halogenated alkanes': [1.2, -0.2],
                           'halogenated alkanes QSPR': [1.22, -0.12],
                           'halogenated alkenes': [1.0 , 0.0],
                           'trihalo-methanes': [1.0, 0.0],
                           'aromatics': [0.9, 0.1],
                           'nitro compounds': [0.75, 0.25],
                           'chlorinated hydrocarbon': [0.59, 0.41],
                           'phenols': [0.65, 0.35],
                           'PNAs': [0.32, 0.68],
                           'pesticides': [0., 0.05],
                           'PFAS': {'PFBA': [0.82, 0.12],
                                    'PFPeA': [0.67, 0.19],
                                    'PFHxA': [0.55, 0.28],
                                    'PFHpA': [0.44, 0.36],
                                    'PFOA': [0.34, 0.44],
                                    'PFNA': [0.24, 0.53],
                                    'PFDA': [0.17, 0.61],
                                    'PFBS': [0.68, 0.24],
                                    'PFHxS': [0.44, 0.48],
                                    'PFHpS': [0.34, 0.57],
                                    'PFOS': [0.25, 0.66],
                                    'PFMOAA': [0.88, 0.055],
                                    'PFO2HxA': [0.63, 0.17],
                                    'PFO3OA': [0.43, 0.27],
                                    'GenX': [0.50, 0.34],
                                    'PFO4DA': [0.25, 0.37],
                                    'NafionBP2': [0.29, 0.63],
                                    '62FTS': [0.38, 0.62],
                                    '82FTS': [0.20, 0.79],
                                    ## Extra, extrapoloated/estimated from above
                                    'PFPrS': [0.63, 0.14], # linear extrapolation for sulfonated C3
                                    'PMPA': [0.63, 0.17],  # copied from PFO2HxA, similar #C #F ether
                                    'PFPeS': [0.49, 0.35], ## linear interpolation for sulfonated C5
                                    'PFPrA': [0.88, 0.03], ## linear interpolation for PFCAs C3
                                    # special added for QSPR Paper
                                    'TCE': [1.0, 0.],
                                    'PFOS_sp': [1.22, -0.12],
                                    'PFOA_sp': [1.22, -0.12],
                                    #average, do not use generally
                                    'Ave': [0.45, 0.41]} 
                           }}
#Update PFHpS

# =============================================================================
# Helper Functions
# =============================================================================

def calc_solver_matrix(nr, nz, ne):
    '''
    nr: number of radial points
    nz: number of axial points
    ne: number of axial segments
    '''
    solver_data = {'nc': nr, 'mc': nz}
    def dspoly(ia, nd, nc, x):
        q = np.ones((nd, nc))
        y = np.zeros(2*nc)
        c = np.zeros((nd, nc))
        d = np.zeros((nd, nc))
        f = np.zeros(nd)
        
        loopq = [2*(i+1)-3 for i in range(1,nc)]
        loopc = [2*(i+1)-4 for i in range(1,nc)]
        cfac = np.array([2*i for i in range(1,nc)])
        loopd = loopq[:-1]
        dfac = np.array([(2*i-2)*(2*i-4+ia) for i in range(3,nc+1)])
    
        for j in range(nd):
            y = np.array([x[j]**n for n in range(1,2*nc+1)])
            q[j,1:] = np.array([y[i] for i in loopq])
            c[j,1:] = cfac * np.array([y[i] for i in loopc])#2. * y[0]
            d[j][1] = 2. * ia
            d[j,2:] = dfac*np.array([y[i] for i in loopd])
            f[j] = 1./float(2.* (j+1) - 2. + ia)
            
        return q, f, c, d 
    
    def dupoly(nd, nc, x):
        q = np.ones((nd, nc))
        y = np.zeros(2*nc)
        c = np.zeros((nd, nc))
        d = np.zeros((nd, nc))
        
        for j in range(nd):
            y[0] = 1.
            for i in range(1,nc):
                y[i] = y[i-1] * (2. * x[j] - 1.)
            q[j,:] = np.array([y[i] for i in range(nc)])
            c[j][1] = 2.
            for i in range(2,nc):
                c[j][i] = 2. * i * y[i-1]
                d[j][i] = 4. * i * (i - 1.) * y[i-2]
        
        return q, c, d
    
    def build_colloc(n_pts):
        '''
        levih method
        '''
        n_pts_m = n_pts - 2 #number of interior points
        
        roots = np.zeros([n_pts])
        roots[-1] = 1 #right boundary
        roots[1:-1] = special.roots_sh_legendre(n_pts_m)[0]
        
        Qmat = np.zeros([n_pts, n_pts])
        Cmat = np.zeros([n_pts, n_pts])
        Dmat = np.zeros([n_pts, n_pts])
        
        for i in range(n_pts):
            for j in range(n_pts):
                Qmat[i, j] = roots[i]**j
            for j in range(1, n_pts):
                Cmat[i, j] = (j) * roots[i]**[j-1]
            for j in range(2, n_pts):
                Dmat[i, j] = (j)*(j-1) * roots[i]**[j-2]
        
        Qinv = np.linalg.inv(Qmat)
        Amat = np.dot(Cmat, Qinv)
        Bmat = np.dot(Dmat, Qinv)
        
        return (roots, Amat, Bmat)
    
    def advect_operator(ne, nz):
        nz_n = ne * nz - (ne - 1) #total number of axial points
        f = ne #Element width adjustment to derivatives
        
        roots, A, _ = build_colloc(nz) # construct first derivative operator
        
        ### Calculation location of overall gridpoints
        Xvals = np.zeros(ne * nz)
        for k in range(ne):
            Xvals[k*nz:(k+1)*nz] = (roots/f + k/f)#/f
        
        to_delete = [nz*(k+1) for k in range(ne-1)]
        Xvals = np.delete(Xvals, to_delete)
        
        Q = A[-1,-1] - A[0,0]
        
        #adjusted element advection operator
        Z = np.zeros([nz, nz])
        Z[:, :] = f * A
        
        #construct overall advection operator
        Adv_Op = np.zeros([nz_n, nz_n])
        
        #fill in blocks for element interiors
        for k in range(ne):
            ii = k*(nz - 1)
            iii = (k+1)*nz - k
            Adv_Op[ii:iii, ii:iii] = Z[:,:]
        
        #fill in continuation points where elements meet
        for k in range(ne-1):
            idx = k * nz - k
            ii = (k + 1)*nz - (k + 1)
            iii = (k + 2)*nz - (k + 1)
        
            Adv_Op[ii, :] = 0 # zero out continuation 
            
            CC1 = Z[-1,:-1] - Z[-1,-1]/Q * A[-1,:-1]
            CC2 = Z[-1, -1]/Q * A[0, 1:]
            
            Adv_Op[ii, idx:ii] = CC1
            Adv_Op[ii, ii+1:iii] = CC2
            
        Adv_Op[0, :] = 0 # Constant inlet boundary
        
        return Adv_Op, Xvals
    
    ngeor = 3 # spherical
    alfar = 1.
    betar = 0.5
    
#   radial matrices
    rr = betar*(special.j_roots(nr-1,alfar,betar)[0]+alfar)
    rr = np.append(rr, [1])#[0],rr)
    r = np.sqrt(rr)
    q, f, c, d = dspoly(ngeor, nr, nr, r)
    qi2 = np.linalg.inv(q)
    wr = f.dot(qi2)#.dot(f)
    br = d.dot(qi2)
    
#    axial
    if ne == 1:
        z = special.roots_sh_legendre(nz-2)[0]
        z = np.append(np.append([0],z),[1])
        q, c, d = dupoly(nz, nz, z)
        qi = np.linalg.inv(q)
        az = c.dot(qi)
    elif ne > 1:
        az, _ = advect_operator(ne, nz)
        solver_data['mc'] = ne * nz - (ne - 1)
    
#save data to export object
    solver_data['wr'] = wr
    solver_data['az'] = az
    solver_data['br'] = br
    
    return solver_data

def ohashi(re, sc):
    ''' returns sh value based on re and sc numbers
    '''
    if re<5.8:
        sh = (2. + 1.58*(re**0.4)*(sc**(1./3.)))
    elif (re>=5.8) and re<500:
        sh = (2. + 1.21*(re**0.5)*(sc**(1./3.)))
    else:
        sh = (2. + 0.59*(re**0.6)*(sc**(1./3.)))
    return sh

def gnielinski(re, sc, ebed):
    ''' returns sh value based on re, sc and ebed numbers'''
    sh_l = 0.644*(re**0.5)*(sc**(1./3.))
    sh_t = (0.037*(re**0.8)*(sc))/(1.+2.443*(re**(-0.1))*((sc**(2./3.))-1.))
    sh = (2. + ((sh_l**2)+(sh_t**2))**(0.5))*(1.+1.5*(1.-ebed))    
    return sh

def fouling_k(k, t, water = 'Rhine', chemical = 'halogenated alkanes'):
    '''
    k : (umole/g)*(L/umole)*(1/n)
    t : minutes
    water= [Rhine, Portage, Karlsruhe, Wausau, Haughton]
    chemical=  [halogenated alkanes, halogenated alkenes, trihalo-methanes
                aromatics, nitro compounds, chlorinated compounds, phenols
                PNAs, pesticides]
    '''
    a1, a2, a3, a4 = foul_params['water'][water]
    b1, b2 = foul_params['chemical'][chemical]
    
    rk1 = b1 * a1 + b2
    rk2 = b1 * a2
    rk3 = b1 * a3
    rk4 = b1 * a4  #no factor of 100, in exponent (should there be b1?)
    
    k_t = k * (rk1 + rk2*t + rk3 * np.exp(rk4 * t))
    
    if type(k_t) == np.ndarray:
        k_t[k_t < k/1e3] = k/1e3
    elif k_t < k/1e3:
        k_t = k/1e3
    return k_t

# tortuosity_slope = 3./(40 * 7 * min_per_day) #40 week interval, 7 days per week.
def tortuosity(t):
    '''
    t in minutes
    
    '''
    if type(t) == np.ndarray:
        tort = 0.334 + 6.61e-6 * t #9.518e-3/min_per_day * t
        tort[tort < 1] = 1.
        return tort
    else:
        if t < 70*min_per_day:
            return 1.
        else:
            return 0.334 + 6.61e-6 * t
        
def viscosity(temp):
    ''' 
    viscosity of water (g/cm-s)
    temp in degC
    
    '''
    t = temp + 273.15
    return np.exp(-24.71 + (4209./t) + 0.04527 * t - (3.376e-5 * t**2))/100.
    
def density(temp):
    '''
    
    density of water gram/cm**3
    temp in degC
    
    '''  
    t = (temp + 273.15)/324.65
    return 0.98396*(-1.41768 + 8.97665*t - 12.2755 * t**2 + 7.45844 * t**3 - 1.73849 * t**4)

def kf_calc(multi_p, re, sc, ebed, corr='Gnielinski'):
#    Working from compiled correlation list from Xu et al. 2013
#       Appl Phys & Eng (2013) 14(3):155-176
    if corr == 'Chern and Chien':
        # Chern and Chien 2002 - modified Gnielinski
        return multi_p*(1.+1.5*(1-ebed))*(2.+0.644*(re**0.5)*(sc**(1./3.))) 
    elif corr == 'Tan et al':
        #Tan et al 1975
        return multi_p*1.1/ebed * (re*sc)**(1./3.)
    elif corr == 'Wilson and Geankoplis':
        #Wilson and Geankoplis, 1966
        return multi_p*1.09/ebed * (re*sc)**(1./3.) 
    elif corr == 'Ohashi et al':
        #Ohashi et al 1981
        return multi_p*ohashi(re, sc), 
    elif corr == 'Williamson et al':
        #Williamson et al 1963
        return multi_p*2.4*ebed*(re**0.3)*(sc**0.42) 
    elif corr == 'Ko et al':
        #Ko et al 2003
        return multi_p*(0.325/(ebed*(re**0.36)*(sc**(1./3.))))
    elif corr == 'Wakau and Funazkri':
        #Wakau and Funazkri, 1978
        return multi_p*(2.+1.1*(re**0.6)*(sc**(1./3.)))
    elif corr == 'Kakaoka et al':
        #Kataoka et al 1972
        return multi_p*(1.85*((1-ebed)/ebed)**(1./3.)*(re**(1./3.))*\
                        (sc**(1./3.)) )
    elif corr == 'Gnielinski':
        #Gnielinski, 1978
        return multi_p*gnielinski(re, sc, ebed) 
    else:
        print('Correlation not found, using Gnielinski')
        return multi_p*gnielinski(re, sc, ebed)


def filter_compounds(df, compounds, carbon, influent):
    #allows dropping of compounds within pandas and drops the column headers
    orig_cols = df.columns.levels[1]
    drop_cols = list(set(orig_cols).difference(compounds))
    tmp_df = df.drop(labels=drop_cols, level=1, axis=1)[[carbon, influent]]
    all_lvls = tmp_df.columns.levels
    codes = tmp_df.columns.codes
    upper = []
    for i in codes[0]:
        upper.append(all_lvls[0][i])
    lower = []
    for i in codes[1]:
        lower.append(all_lvls[1][i])
    new_lvls = [upper] + [lower]
    
    m_idx = pd.MultiIndex.from_arrays(new_lvls, names=['type','compound'])
    
    tmp_df.columns=m_idx
    return tmp_df

def recalc_k(orig_k, q, xn_old, xn_new):
    return q**(1-(xn_new/xn_old))*orig_k**(xn_new/xn_old)

def interp(array, delta):
    #assumes delta between array of 1
    #array is a [0, 1] dimensional array
    return (array[1]-array[0]) * delta + array[0]

def calc_Jac(ff, u0, hh):
    """
    FUNCTION IS SLOW
    USE: spar_Jac instead. It DOES NOT require running diffun.
    
    calculate jacobian of ff on u using finite differences with step size hh
    
    Similar to how solve_ivp solves the same problem. 
    Requires running the diffun "len(u0)" times
    
    """
    NEQ = len(u0)
    Jac = np.zeros((NEQ, NEQ))
    ff0 = ff(0, u0)

    for ii in range(NEQ):  # XXX: slow and with lots of array copies. 
        u_high = u0 + 0.0
        u_high[ii] = u0[ii] + hh
        Jac[ii, :] = (ff(0, u_high) - ff0)/(hh)
        
    return Jac.T

def spar_Jac(ncomp, nr, nz, ne):
    '''
    
    Helper function that predetermines the sparsity matrix. This is provided to 
    solve_ivp to prevent solve_ivp from having to solve for this itself. solve_ivp
    runs the function's diffun to solve this problem, which is much slower. 
    
    Called in PSDM.__init__()
    
    Parameters
    ----------
    ncomp : int
        number of compounds to consider.
    nr : int
        number or radial collocation points.
    nz : int
        number or axial collocation points.
    ne : int
        number of axial elements.

    Returns
    -------
    sparse matrix that indicates sparsity to solve_ivp.

    '''
    nz_tot = (nz * ne) - (ne - 1)    
    cube_size = nz_tot * (nr + 1) 
    bead_num = nr * nz_tot
    
    #Consistent features
    diags = np.ones((2*(nr-1)+ 1, bead_num))
    diag_loc=np.arange(-((nr-1)*nz_tot), (nr-1)*nz_tot+1, nz_tot)

    #bulk diagonal matrix
    mat1 = sparse.spdiags(diags, diag_loc, bead_num, bead_num)
    
    #always true
    low_diag = np.ones(bead_num) #might be nz_tot-1
    low_diag[-nz_tot] = 0. #drops out if not nz_tot
    lower_mat = sparse.spdiags(low_diag, (nr-1)*(nz_tot), nz_tot, bead_num)
       
    mat2 = sparse.bmat([[mat1],[lower_mat]])
    right1 = sparse.bsr_matrix((bead_num-nz_tot, nz_tot))

    #features that change with sub-location
    for i in range(ncomp):
        for j in range(ncomp):
            if i==j:
                rightp = np.ones(nz_tot)
                rightp[0] = 0.
                right2 = sparse.diags(rightp, 0)
                right3 = np.ones((nz_tot, nz_tot))
                right3[0] = 0.
                right3[:,0] = 0.
            else:
                right2 = sparse.bsr_matrix((nz_tot, nz_tot))
                right3 = right2
            right_mat = sparse.bmat([[right1],[right2]])
            right_mat = sparse.bmat([[right_mat], [right3]])
            
            mat1 = sparse.bmat([[mat2, right_mat]])
        
            if j==0:
                rowmat = mat1
            else:
                rowmat = sparse.bmat([[rowmat, mat1]])
        if i==0:
            newmat = rowmat
        else:
            newmat = sparse.bmat([[newmat],[rowmat]])
            
    
    return newmat # mat1, mat2
        
def find_minimum_df(df):
    def get_min(df):
        '''
        convenience function
        '''
        return (df[df==df.min().min()]).dropna(how='all').dropna(how='all', axis=1)
    
    if isinstance(df, pd.DataFrame):
        max_val = df.max().max()
        thresh_val = np.percentile(df,10)
                
        np_data = df.values * 1.
        numX, numY = np.shape(np_data)
        
        oldX, oldY = np.where(np_data == max_val)
        
        val_list = []
        
        orig_min_val = df.min().min()
        orig_minX, orig_minY = np.where(df.values == orig_min_val)
        orig_min = pd.DataFrame(orig_min_val, \
                                columns = [df.columns[int(orig_minY)]],\
                                index = [df.index[int(orig_minX)]])
        
        for i in range(10):#while not done:
            stepsX = 3
            stepsY = 3
            
            cur_min = np.min(np_data)
            xloc, yloc = np.where(np_data == cur_min)
            
            if xloc[0] > 0:
                minX = xloc[0] - 1
            else:
                minX = 0
            
            if xloc[0] < numX - 1:
                maxX = xloc[0] + 2 # ?
            else:
                maxX = numX
                
            if yloc[0] > 0:
                minY = yloc[0] - 1
            else:
                minY = 0
            
            if yloc[0] < numY - 1:
                maxY = yloc[0] + 2 # ?
            else:
                maxY = numY
            
            if xloc[0] == 0 or xloc[0] == numX-1:
                stepsX = 2
            if yloc[0] == 0 or yloc[0] == numY-1:
                stepsY = 2
            
            tol_level = stepsX * stepsY #might consider -1
            subset = np_data[minX:maxX,minY:maxY]
            
            test_mean = np.mean(subset - cur_min)/cur_min
            
            val_list.append([test_mean, cur_min, \
                             xloc[0], yloc[0], \
                             np.sum(subset<thresh_val), tol_level])

            np_data[xloc, yloc] = max_val
            
        np_val = np.array(val_list)
        best_min = np.min(np_val[:,0])
        best_loc = np.where(best_min == np_val[:,0])
        bestM, bestMin, bestX, bestY, bestThresh, actTol = np_val[best_loc[0]][0]
        
        idx = [df.index[int(bestX)]]
        col = [df.columns[int(bestY)]]
        
        if bestThresh == actTol:
            df_out = pd.DataFrame(bestMin, columns=col, index=idx)
        else:
            df_out = orig_min.copy()
        
    else:
        df_out = pd.DataFrame(0, columns=[1],index=1)
        print('Error: find_minumum_df must be provided a dataframe')
    
    if df.values.shape != (1,1):
        if len(df_out.columns) != 1 or len(df_out.index) != 1:
            df_out = orig_min.copy()
        
    return df_out

def generate_grid(grid_num_xn, grid_num_k,\
                  loop_num=0, xn_range=[1], k_range=[1]):
    
    #generates optimizer grid for staged approach
    #called by run_all
    if loop_num == 0 and len(xn_range) == 1:
        max_xn = 0.9
        min_xn = 0.3
    else:
        max_xn = np.max(xn_range)
        min_xn = np.min(xn_range)
    
    xn_range_new = np.linspace(min_xn, max_xn, grid_num_xn)
    # grid_xn = (max_xn - min_xn)/(grid_num_xn - 1)

    if loop_num == 0 and len(k_range)==1:
        max_k = 5.
        min_k = 1.
    else:
        max_k = np.max(k_range)
        min_k = np.min(k_range)
    k_range_new = np.linspace(min_k, max_k, grid_num_k)
    # grid_k = (max_k - min_k)/(grid_num_k - 1)    

    return xn_range_new, k_range_new




def logistic(t, c0, a, b):
    ## Returns Logistic form function
    return c0 / (1 + np.exp(a - b * t))









