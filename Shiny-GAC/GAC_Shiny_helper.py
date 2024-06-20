
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

import sys
sys.path.append('../../Water_Treatment_Models/PSDM/')


import PSDM





def run_PSDM(columndata, chem_data, kdata, infdat, effdat, nr, nz, water_type='Organic Free', chem_type='halogenated alkenes'): ## need fouling passed from R as well
    ## convert from R to Python need
    column_info = pd.Series(columndata['value'].values, index=columndata['name'])
    
   
    
    for i in ['rad', 'epor', 'psdfr', 'rhop', 'rhof', 
              'L', 'wt', 'flrt', 'diam', 'tortu', 
              'mass_mul', 'flow_mult', 't_mult']: # changing numeric values to floats
        column_info[i] = pd.to_numeric(column_info[i]).astype(float) 
       
    ## PREPARE chemical property data
    prop_columns = list(chem_data.columns) ### compounds
    compounds = prop_columns[1:]
    chem_updated = pd.DataFrame(chem_data[compounds].values, #.astype(float),
                                columns=compounds,
                                index=chem_data[prop_columns[0]].values)
    
    # ## PREPARE k_data information
    # ## assumes shape of the data is the same
    kdata_updated = pd.DataFrame(kdata[compounds].values,
                                 columns=compounds,
                                 index=kdata[prop_columns[0]])

    midx = [(i, j) for i in [column_info['influentID'], column_info['effluentID']] for j in prop_columns[1:]]

    max_inf_time = np.max(infdat['time'].values) ### only consider data in influent
    interval = 1

    raw_data = pd.DataFrame(0, columns=pd.MultiIndex.from_tuples(midx), index=np.arange(0, max_inf_time+interval, interval))

    for comp in compounds:
        f_infl = interp1d(infdat['time'].values, infdat[comp].values)
        f_effl = interp1d(effdat['time'].values, effdat[comp].values, fill_value='extrapolate')
        for time in raw_data.index:
            raw_data.loc[time][column_info['influentID'], comp] = f_infl(time)
            raw_data.loc[time][column_info['effluentID'], comp] = f_effl(time)

    
 
    
    ## set up output dataframe for consistency
    output_data = pd.DataFrame(0, columns=compounds, index=np.arange(0,max_inf_time+.2, 0.25))

    try:
        column = PSDM.PSDM(column_info,
                               chem_updated,
                               raw_data,
                               k_data=kdata_updated,
                               nr=int(nr),
                               nz=int(nz),
                               water_type=water_type,
                               chem_type=chem_type
                               )
                               

        
        column.model_uncertainty(capacity='None') ## uses model_uncertainty to handle running sim
        
        for comp in compounds:
            f_compound = interp1d(column.results.index, 
                                  column.results[comp].values,
                                  fill_value='extrapolate')
            output_data[comp] = f_compound(output_data.index)

    except Exception as e:
        return e
      
    output_data['time'] = output_data.index
    return output_data
  
  #### Function for integrating fitting capability
def run_PSDM_fitter(columndata, chem_data, kdata, infdat, effdat, nr, nz, water_type, chem_type, 
                      pm=30, xn=0.01): ## need fouling passed from R as well #pm=30, xn=0.01
    ## convert from R to Python need
    column_info = pd.Series(columndata['value'].values, index=columndata['name'])
   
    
    for i in ['rad', 'epor', 'psdfr', 'rhop', 'rhof', 
              'L', 'wt', 'flrt', 'diam', 'tortu', 
              'mass_mul', 'flow_mult', 't_mult']: # changing numeric values to floats
        column_info[i] = pd.to_numeric(column_info[i]).astype(float) 

       
    ## PREPARE chemical property data
    prop_columns = list(chem_data.columns) ### compounds
    compounds = prop_columns[1:]
    chem_updated = pd.DataFrame(chem_data[compounds].values, #.astype(float),
                                columns=compounds,
                                index=chem_data[prop_columns[0]].values)
    
    
    # ## PREPARE k_data information
    # ## assumes shape of the data is the same
    kdata_updated = pd.DataFrame(kdata[compounds].values,
                                 columns=compounds,
                                 index=kdata[prop_columns[0]])

    midx = [(i, j) for i in [column_info['influentID'], column_info['effluentID']] for j in prop_columns[1:]]

    max_inf_time = np.max(infdat['time'].values) ### only consider data in influent
    interval = 1

    raw_data = pd.DataFrame(0, columns=pd.MultiIndex.from_tuples(midx), index=np.arange(0, max_inf_time+interval, interval))

    for comp in compounds:
        f_infl = interp1d(infdat['time'].values, infdat[comp].values)
        f_effl = interp1d(effdat['time'].values, effdat[comp].values, fill_value='extrapolate')
        for time in raw_data.index:
            raw_data.loc[time][column_info['influentID'], comp] = f_infl(time)
            raw_data.loc[time][column_info['effluentID'], comp] = f_effl(time)

    
    ## set up output dataframe for consistency
    output_data = pd.DataFrame(0, columns=compounds, index=np.arange(0,max_inf_time+.2, 0.25))

    try:
        column = PSDM.PSDM(column_info,
                               chem_updated,
                               raw_data,
                               k_data=kdata_updated,
                               nr=int(nr),
                               nz=int(nz),
                               water_type=water_type,
                               chem_type=chem_type
                               )
                               

        column.run_all_smart(pm=pm, save_file=False, num=11)
        
        column.model_uncertainty(capacity='None') ## uses model_uncertainty to handle running sim
        
        for comp in compounds:
            f_compound = interp1d(column.results.index, 
                                  column.results[comp].values,
                                  fill_value='extrapolate')
            output_data[comp] = f_compound(output_data.index)

    except Exception as e:
        return e, []
      
    output_data['time'] = output_data.index
    return output_data, column.k_data  ## returns model fits and associated k_data table.
