# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 08:48:47 2021

NOTE: This test is only valid for monovalent trace contaminants.

Also note that this is in the film limited regime.
We might want to add a check for HSDM Deff for PSDM model. 

@author: LHaupert
"""

from context import hsdmix, psdmix

import os
import unittest

import numpy as np
from scipy.special import erf

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


def HSDM_ana(T, Dg, Bi, St):
    """
    Rosen's asymptotic solution for linear isotherm.
    Should work for a trace contaminant in ion exchange.
    
    Returns C/C0
    """
    
    term1 = ((Dg + 1) / Dg) * (T - 1)  # arguably easier to understand
    term2 = 2 * np.sqrt((5 + Bi) / (15 * St)) 
    C_C0 = 0.5 * (1 + erf(term1 / term2))
    
    return C_C0



hours = np.linspace(0, 200, num=100)

HSDM = hsdmix.HSDMIX(DATA_DIR+'/trace_test_hsdmix.xlsx')

vs = HSDM.params['v']
L = HSDM.params['L']
Ds = HSDM.params['Ds']
kL = HSDM.params['kL']
EBED = HSDM.params['EBED']
rb = HSDM.params['rb']

Q = HSDM.params['Q']
alphaBA = HSDM.ions.iloc[1, 1]
CA = HSDM.Cin_t.iloc[0, 0]
CB = HSDM.Cin_t.iloc[0, 1]

lambda_B = Q / (CB + CA / alphaBA)

vi = vs / EBED # interstitial velocity (cm/s)
tau = L / vi # Packed bed contact time
Dg = lambda_B * (1 - EBED) / EBED # mass distribution coeff (HSDM)
Bi = (kL * rb * (1 - EBED) / (Ds * Dg * EBED)) # Biot number
St = kL * tau * (1 - EBED) / (EBED * rb) # Stanton number
Ed = Ds * Dg * tau / (rb**2) # diffusion modulus

t_eval = hours * 60 * 60 # (seconds)
T_eval = t_eval / (tau * (Dg + 1))


y_Rosen = HSDM_ana(T_eval, Dg, Bi, St)


class TraceTestSuite(unittest.TestCase):
    """Test HSDMIX and PSDMIX models for trace contaminant."""
    

    def test_trace_HSDMIX(self):
        """ Simple trace contaminant (HSDMIX) """
        IEX = hsdmix.HSDMIX(DATA_DIR+'/trace_test_hsdmix.xlsx')
        t, u = IEX.solve(t_eval=hours, const_Cin=True)
 
        assert np.allclose(u[0, -1, -1, :]/CB, y_Rosen, atol=0.03)

    def test_trace_PSDMIX(self):
        """ Simple trace contaminant (PSDMIX) """
        IEX = hsdmix.HSDMIX(DATA_DIR+'/trace_test_psdmix.xlsx')
        t, u = IEX.solve(t_eval=hours, const_Cin=True)
 
        assert np.allclose(u[0, -1, -1, :]/CB, y_Rosen, atol=0.03)
   

if __name__ == '__main__':
    unittest.main()






