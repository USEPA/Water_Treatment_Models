# -*- coding: utf-8 -*-

from context import psdmix

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
u_check = np.array([6.99991177, 6.46910715, 5.11604299, 3.92926336, 3.11649178,
                    2.60910929, 2.3143242 , 2.15397154, 2.06881489, 2.01618331,
                    1.96902224, 1.91151521, 1.83819327, 1.75624193, 1.68180112,
                    1.62738689, 1.59348207, 1.57354937, 1.56106999, 1.55137199,
                    1.54265415, 1.53466305, 1.52789321, 1.52255164, 1.51852714,
                    1.51562096, 1.51363523, 1.51232997, 1.51148386, 1.51093083,
                    1.51057451, 1.51034814, 1.51022782, 1.51015245, 1.51011182,
                    1.51008239, 1.51005762, 1.51004523, 1.51002954, 1.51001285,
                    1.50999839, 1.50999199, 1.50998687, 1.50998303, 1.50998046,
                    1.50997918, 1.50997917, 1.50998044, 1.50998324])


class PSDMIXTestSuite(unittest.TestCase):
    """Test PSDMIX model."""
    
    def setUp(self):
        warnings.simplefilter('ignore', PendingDeprecationWarning)

    def test_regress(self):
        """ Simple Regression Test """
        IEX = psdmix.PSDMIX(DATA_DIR+'/reg_test_input_PSDM.xlsx')
        t, u = IEX.solve(t_eval=hours, const_Cin=True)
 
        assert np.allclose(u[0,0,-1,:], u_check, rtol=0.01)


   

if __name__ == '__main__':
    unittest.main()
