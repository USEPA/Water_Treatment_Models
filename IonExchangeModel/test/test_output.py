# -*- coding: utf-8 -*-

from context import hsdmix

import io
import os
import unittest
import warnings

import numpy as np
import pandas as pd


class Results(object):
    def __init__(self, t, y):
        self.t = t
        self.y = y

t = np.zeros(2)
y = np.zeros((16, 2))
result = Results(t, y)

u_result = np.zeros((2, 4, 2, 2))
u_result[0, :, -1, :] = 1.51


DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

class OutputTestSuite(unittest.TestCase):
    """Test Collocation construction."""

    def setUp(self):
        warnings.simplefilter('ignore', PendingDeprecationWarning)
        self.infile = DATA_DIR+'/reg_test_input.xlsx'


    def test_out_units_mg(self):

        ions = pd.read_excel(self.infile, sheet_name='ions', index=0)
        Cin = pd.read_excel(self.infile, sheet_name='Cin', index=0)        
        params = pd.read_excel(self.infile, sheet_name='params', index=0)
        
        ions.loc[0, 'units'] = 'mg'
        Cin.iloc[:, 1] = Cin.iloc[:, 1] * ions.loc[0, 'mw']
        
        mock_infile = io.BytesIO()
        writer = pd.ExcelWriter(mock_infile, engine='xlsxwriter')
        params.to_excel(writer, sheet_name='params', index=0)
        ions.to_excel(writer, sheet_name='ions', index=0)
        Cin.to_excel(writer, sheet_name='Cin', index=0)    
        writer.save()
        
        IEX = hsdmix.HSDMIX(writer) 
       
        writer.close()
        
        # To save time, assign IEX fields using prebuilt arrays 
        # instead of solving.
        IEX.result = result
        IEX.timeback = 1
        IEX.u_result = u_result
        
        mock_outfile = io.BytesIO()
        IEX.save_results(mock_outfile)
        
        result_df = pd.read_excel(mock_outfile, index=0)
        Cl_out = result_df.values[-1, 1]
                
        assert np.allclose(Cl_out, 53.52952053)
        

if __name__ == '__main__':
    unittest.main()
