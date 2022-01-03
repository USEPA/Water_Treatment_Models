# -*- coding: utf-8 -*-

from context import hsdmix

import io
import os
import unittest
import warnings

import numpy as np
import pandas as pd


DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

class InputTestSuite(unittest.TestCase):
    """Test Collocation construction."""

    def setUp(self):
        warnings.simplefilter('ignore', PendingDeprecationWarning)
        self.infile = DATA_DIR+'/reg_test_input.xlsx'


    def test_in_units_mg(self):

        ions = pd.read_excel(self.infile, sheet_name='ions', index_col=0)
        Cin = pd.read_excel(self.infile, sheet_name='Cin', index_col=0)        
        params = pd.read_excel(self.infile, sheet_name='params', index_col=0)
        
        ions.loc['CHLORIDE', 'units'] = 'mg'
        Cin.iloc[:, 0] = Cin.iloc[:, 0] * ions.loc['CHLORIDE', 'mw']
        
        mock_infile = io.BytesIO()
        writer = pd.ExcelWriter(mock_infile, engine='xlsxwriter')
        params.to_excel(writer, sheet_name='params')
        ions.to_excel(writer, sheet_name='ions')
        Cin.to_excel(writer, sheet_name='Cin')    
        writer.save()
        
        IEX = hsdmix.HSDMIX(writer.handles.handle) 
       
        writer.close()
                
        assert np.allclose(IEX.Cin_t.iloc[:, 0], np.array([1.51, 1.51]))
        

if __name__ == '__main__':
    unittest.main()
