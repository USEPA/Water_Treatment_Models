# -*- coding: utf-8 -*-

from context import hsdmix

import os
import unittest
import warnings

import numpy as np


DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

hours = np.linspace(0, 48, num=49)
# Check the presaturant at the outlet.
# Avoids having to match numerical noise at the start of run.
# Presaturant outlet is usually sensitive to problems in other areas.
# However, more rigorous testing may be required eventually.  
u_check = np.array([7.        , 5.98772214, 4.81410572, 3.91171217, 3.25520516,
                    2.79727629, 2.48876065, 2.28671097, 2.15422175, 2.06077632,
                    1.98178017, 1.8993725 , 1.80707225, 1.71343639, 1.63558986,
                    1.58273139, 1.55116331, 1.53341656, 1.52364305, 1.51847861,
                    1.51546164, 1.5135453 , 1.51234699, 1.51160103, 1.51111173,
                    1.51075915, 1.5105166 , 1.51035164, 1.5102427 , 1.51017731,
                    1.51012353, 1.51009665, 1.51007347, 1.51005018, 1.51003757,
                    1.51002372, 1.51001366, 1.51000719, 1.51000234, 1.50999681,
                    1.50999288, 1.50999054, 1.50998979, 1.50999063, 1.50999307,
                    1.50999374, 1.50999404, 1.50999466, 1.50999561])


class HSDMIXTestSuite(unittest.TestCase):
    """Test HSDMIX model."""
    
    def setUp(self):
        warnings.simplefilter('ignore', PendingDeprecationWarning)

    def test_regress(self):
        """ Simple Regression Test """
        IEX = hsdmix.HSDMIX(DATA_DIR+'/reg_test_input.xlsx')
        t, u = IEX.solve(t_eval=hours, const_Cin=True)
 
        assert np.allclose(u[0,0,-1,:], u_check, rtol=0.01)

    def test_regress_alt(self):
        """ Alternate Regression Test """
        IEX = hsdmix.HSDMIX(DATA_DIR+'/reg_test_input_alt.xlsx')
        t, u = IEX.solve(t_eval=hours, const_Cin=True)
 
        assert np.allclose(u[0,0,-1,:], u_check, rtol=0.01)

   

if __name__ == '__main__':
    unittest.main()
