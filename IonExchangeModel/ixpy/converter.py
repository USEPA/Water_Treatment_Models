import pandas as pd
import numpy as np

min_per_day = 24 * 60
s_per_day = min_per_day * 60

def conv_iex_u(u, t, ions):
    ## Converter function to produce easier to handle DataFrames
    # converts from standard model output structure, returns time indexed DataFrame, 
    # converted to input units in concentration by column

    ix_df = pd.DataFrame(index=np.round(t/s_per_day, 2), ## rounding prevents indexing error in loop
                         columns=ions.index)

    num_compounds = len(ions.index)
    for i in range(num_compounds):
        factor = 1 ## assumes meq 
        comp_name = ions.index[i]

        if ions.loc[comp_name]['units'] == 'ug':
            factor = ions.loc[comp_name]['mw'] * 1e3 ## meq -> ug/L
        elif ions.loc[comp_name]['units'] == 'ng':
            factor = ions.loc[comp_name]['mw'] * 1e6 ## meq -> ng/L
        elif ions.loc[comp_name]['units'] == 'mg':
            factor = ions.loc[comp_name]['mw'] * 1 ## meq -> mg/L
        else: 
            factor = 1000 / (ions.loc[comp_name]['mw'] * ions.loc[comp_name]['valence']) 
        ### what about mgC, mgN units??? 
        ## TODO: >>
        
        ix_df[comp_name] = u[0, i, -1, :] * factor ### returned as meq... need to add smarter conversion
        
        
        # else:
            # ix_df[ions['name'][i]] = u[0, i, -1, :]
        
    ix_df[ix_df <= 0] = 1e-16 ## make negatives or 0's very small, numerical stability

    return ix_df
