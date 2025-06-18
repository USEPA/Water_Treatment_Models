import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

import sys
sys.path.append('./psdm/')

import PSDM

python = False ## to all testing in Jupyter Notebook rather than in R for debugging issues

def run_PSDM(columndata, chem_data, kdata, infdat, effdat, nr, nz, water_type='Organic Free', chem_type='halogenated alkenes'): ## need fouling passed from R as well
    ## convert from R to Python need
    column_info = pd.Series(columndata['value'].values, index=columndata['name'])

    ## PREPARE chemical property data
    prop_columns = list(chem_data.columns) ### compounds
    #print(prop_columns)
    mass_transfer_df = 'None'
    if python: ## testing python
        compounds = prop_columns[0:]
        chem_updated = chem_data.copy()

        # ## PREPARE k_data information
        # ## assumes shape of the data is the same
        kdata_updated = kdata.copy()

    else:
        for i in ['rad', 'epor', 'psdfr', 'rhop', 'rhof',
              'L', 'wt', 'flrt', 'diam', 'tortu',
              'mass_mul', 'flow_mult', 't_mult']: # changing numeric values to floats
            column_info[i] = pd.to_numeric(column_info[i])


        compounds = prop_columns[1:]
        chem_updated = pd.DataFrame(chem_data[compounds].values, #.astype(float),
                                    columns=compounds,
                                    index=chem_data[prop_columns[0]].values)
        
        mass_transfer_df = pd.DataFrame(0, index=['kf', 'dp', 'ds'], columns=compounds)

        for idx in chem_updated.index:
            try:
                lower_idx = idx.lower()
                for comp in compounds:
                    if lower_idx in ['kf', 'dp', 'ds']:
                        mass_transfer_df.loc[lower_idx, comp] = chem_updated.loc[idx, comp] * 1
            except:
                pass


        # ## PREPARE k_data information
        # ## assumes shape of the data is the same
        kdata_updated = pd.DataFrame(kdata[compounds].values,
                                     columns=compounds,
                                     index=kdata[prop_columns[0]])

    midx = [(i, j) for i in [column_info['influentID'], column_info['effluentID']] for j in compounds]

    max_inf_time = np.max(infdat['time'].values) ### only consider data in influent
    interval = 1

    raw_data = pd.DataFrame(0, columns=pd.MultiIndex.from_tuples(midx), index=np.arange(0, max_inf_time+interval, interval))

    idx = pd.IndexSlice
    for comp in compounds:
        raw_data.loc[:, idx[column_info['influentID'], comp]] = np.interp(raw_data.index, infdat['time'].values, infdat[comp].values)
        raw_data.loc[:, idx[column_info['effluentID'], comp]] = np.interp(raw_data.index, effdat['time'].values, effdat[comp].values)

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
                               chem_type=chem_type,
                               mass_transfer=mass_transfer_df
                               )

        column.model_uncertainty(capacity='None') ## uses model_uncertainty to handle running sim

        for comp in compounds:
            output_data[comp] = np.interp(output_data.index,
                                          column.results.index,
                                          column.results[comp].values)

    except Exception as e:
        return e

    output_data['time'] = output_data.index
    return output_data

  #### Function for integrating fitting capability
def run_PSDM_fitter(columndata, chem_data, kdata, infdat, effdat, nr, nz, water_type='Organic Free', chem_type='halogenated alkenes', pm=30, xn=0.01): ## need fouling passed from R as well
    ## convert from R to Python need
    column_info = pd.Series(columndata['value'].values, index=columndata['name'])

    ## PREPARE chemical property data
    prop_columns = list(chem_data.columns) ### compounds

    if python: ## testing python
        compounds = prop_columns[0:]
        chem_updated = chem_data.copy()

    else:
        for i in ['rad', 'epor', 'psdfr', 'rhop', 'rhof',
              'L', 'wt', 'flrt', 'diam', 'tortu',
              'mass_mul', 'flow_mult', 't_mult']: # changing numeric values to floats
            column_info[i] = pd.to_numeric(column_info[i])


        compounds = prop_columns[1:]
        chem_updated = pd.DataFrame(chem_data[compounds].values, 
                                    columns=compounds,
                                    index=chem_data[prop_columns[0]].values)


    midx = [(i, j) for i in [column_info['influentID'], column_info['effluentID']] for j in compounds]

    max_inf_time = np.max(infdat['time'].values) ### only consider data in influent
    interval = 1

    raw_data = pd.DataFrame(0, columns=pd.MultiIndex.from_tuples(midx), index=np.arange(0, max_inf_time+interval, interval))

    idx = pd.IndexSlice
    for comp in compounds:
        raw_data.loc[:, idx[column_info['influentID'], comp]] = np.interp(raw_data.index, infdat['time'].values, infdat[comp].values)
        raw_data.loc[:, idx[column_info['effluentID'], comp]] = np.interp(raw_data.index, effdat['time'].values, effdat[comp].values)

    ## set up output dataframe for consistency
    output_data = pd.DataFrame(0, columns=compounds, index=np.arange(0,max_inf_time+.2, 0.25))

    try:
        column = PSDM.PSDM(column_info,
                               chem_updated,
                               raw_data,
                               nr=int(nr),
                               nz=int(nz),
                               water_type=water_type,
                               chem_type=chem_type,
                               )

        try:
            #print(column.k_data)

            column.run_all_smart(pm=float(pm),
                                 des_xn=float(xn),
                                 save_file=False,
                                 num=11)
        except Exception as e:
            returned_data = pd.DataFrame(np.ones(5), index=['K', '1/n', 'q', 'brk','AveC'])
            #print(e, 'Error on run_all_smart')
            return e, returned_data

        column.model_uncertainty(capacity='None') ## uses model_uncertainty to handle running sim

        for comp in compounds:
            # output_data[comp] = np.ones(len(output_data.index))
            output_data[comp] = np.interp(output_data.index,
                                          column.results.index,
                                          column.results[comp].values)

    except Exception as e:
        returned_data = pd.DataFrame(np.ones(5) * 2, index=['K', '1/n', 'q', 'brk','AveC'])
        #print(e, 'Error on column creation')
        return e, returned_data

    output_data['time'] = output_data.index
    return output_data, column.k_data  ## returns model fits and associated k_data table.

