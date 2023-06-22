# -*- coding: utf-8 -*-

from context import hsdmix
from context import psdmix

import os
import unittest
import warnings

import numpy as np


DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

# =============================================================
# Referent values at pH in the range from 8.3 to 6.0
# 100 mg CaCO3 as mg HCO3-  # XXX: meq HCO3- ?
# 2.5 meq CaCO3 as meq HCO3-
# =============================================================

ref_mg = [1.95788, 1.94614, 1.92026, 1.87435, 1.79931,
          1.68378, 1.51777, 1.29989, 1.04502, 0.78395, 0.55042]

ref_meq = [2.44955, 2.43486, 2.40248, 2.34504, 2.25116,
           2.10662, 1.89892, 1.62632, 1.30745, 0.98082, 0.68865]


class AlkalinityTestSuite(unittest.TestCase):
    """Test Total Alkalinity to Bicarbonate ion conversion."""
    
    def setUp(self):
        warnings.simplefilter('ignore', PendingDeprecationWarning)
        self.infile1 = DATA_DIR+'/alkalinity_test_mg.xlsx'
        self.infile2 = DATA_DIR+'/alkalinity_test_meq.xlsx'

    def test_mg_HSDM(self):
        """ HSDM, Total Alkalinity in mg/L as CaCO3 """
        IEX = hsdmix.HSDMIX(self.infile1)
 
        assert np.allclose(IEX.Cin_t['BICARBONATE'], ref_mg)

    def test_meq_HSDM(self):
        """ HSDM, Total Alkalinity in meq/L as CaCO3 """
        IEX = hsdmix.HSDMIX(self.infile2)
 
        assert np.allclose(IEX.Cin_t['BICARBONATE'], ref_meq)

    def test_mg_PSDM(self):
        """ HSDM, Total Alkalinity in mg/L as CaCO3 """
        IEX = psdmix.PSDMIX(self.infile1)
 
        assert np.allclose(IEX.Cin_t['BICARBONATE'], ref_mg)

    def test_meq_PSDM(self):
        """ HSDM, Total Alkalinity in meq/L as CaCO3 """
        IEX = psdmix.PSDMIX(self.infile2)
 
        assert np.allclose(IEX.Cin_t['BICARBONATE'], ref_meq)
   

if __name__ == '__main__':
    unittest.main()
