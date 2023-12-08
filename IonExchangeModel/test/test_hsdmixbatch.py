# -*- coding: utf-8 -*-

from context import hsdmixbatch

import os
import unittest
import warnings

import numpy as np
from numpy import sqrt


DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


def helfferich_film(t, Q, KAB, rb, kL, cAi, cBi, Vw, Vx):
    """
    Return C/C0 over time
    A is contaminant
    B is presaturant
    NOTE: backwards notationally to what is used in our papers
    """

    result = binary_IEX_eq(cAi, cBi, 0, Q, Vw, Vx, KAB)   
    cAf, cBf, qAf, qBf = result 
    
    kAxw = qAf / cAf    # Distribution coefficient
    LAM = kAxw * Vx / Vw
    ex_fact = -3 * kL / rb / kAxw   # eq 30, 31, except for t
    f = 1 - np.exp((LAM+1)*(ex_fact * t)) # eq 67
    C_C0 = 1 - f*(cAi - cAf)/cAi

    return C_C0

def binary_IEX_eq(cAi, cBi, qAi, qBi, Vw, Vx, KAB):
    """
    Solve equilibrium concentrations for batch IX system with
    two monovalent ions.
    """
    Q = qAi + qBi    

    ######## solution 1 ##########
    
    cAf1 = (-KAB*Q*Vx + KAB*Vw*cAi + KAB*Vx*qAi + Q*Vx - 2*Vw*cAi - Vw*cBi - 2*Vx*qAi - Vx*qBi - sqrt(KAB**2*Q**2*Vx**2 - 2*KAB**2*Q*Vw*Vx*cAi - 2*KAB**2*Q*Vx**2*qAi + KAB**2*Vw**2*cAi**2 + 2*KAB**2*Vw*Vx*cAi*qAi + KAB**2*Vx**2*qAi**2 - 2*KAB*Q**2*Vx**2 + 2*KAB*Q*Vw*Vx*cAi + 2*KAB*Q*Vw*Vx*cBi + 2*KAB*Q*Vx**2*qAi + 2*KAB*Q*Vx**2*qBi + 2*KAB*Vw**2*cAi*cBi + 2*KAB*Vw*Vx*cAi*qBi + 2*KAB*Vw*Vx*cBi*qAi + 2*KAB*Vx**2*qAi*qBi + Q**2*Vx**2 - 2*Q*Vw*Vx*cBi - 2*Q*Vx**2*qBi + Vw**2*cBi**2 + 2*Vw*Vx*cBi*qBi + Vx**2*qBi**2))/(2*Vw*(KAB - 1))
    qAf1 = (KAB*Q*Vx + KAB*Vw*cAi + KAB*Vx*qAi - Q*Vx + Vw*cBi + Vx*qBi + sqrt(KAB**2*Q**2*Vx**2 - 2*KAB**2*Q*Vw*Vx*cAi - 2*KAB**2*Q*Vx**2*qAi + KAB**2*Vw**2*cAi**2 + 2*KAB**2*Vw*Vx*cAi*qAi + KAB**2*Vx**2*qAi**2 - 2*KAB*Q**2*Vx**2 + 2*KAB*Q*Vw*Vx*cAi + 2*KAB*Q*Vw*Vx*cBi + 2*KAB*Q*Vx**2*qAi + 2*KAB*Q*Vx**2*qBi + 2*KAB*Vw**2*cAi*cBi + 2*KAB*Vw*Vx*cAi*qBi + 2*KAB*Vw*Vx*cBi*qAi + 2*KAB*Vx**2*qAi*qBi + Q**2*Vx**2 - 2*Q*Vw*Vx*cBi - 2*Q*Vx**2*qBi + Vw**2*cBi**2 + 2*Vw*Vx*cBi*qBi + Vx**2*qBi**2))/(2*Vx*(KAB - 1))
    cBf1 = (-KAB*Q*Vx + KAB*Vw*cAi + 2*KAB*Vw*cBi + KAB*Vx*qAi + 2*KAB*Vx*qBi + Q*Vx - Vw*cBi - Vx*qBi + sqrt(KAB**2*Q**2*Vx**2 - 2*KAB**2*Q*Vw*Vx*cAi - 2*KAB**2*Q*Vx**2*qAi + KAB**2*Vw**2*cAi**2 + 2*KAB**2*Vw*Vx*cAi*qAi + KAB**2*Vx**2*qAi**2 - 2*KAB*Q**2*Vx**2 + 2*KAB*Q*Vw*Vx*cAi + 2*KAB*Q*Vw*Vx*cBi + 2*KAB*Q*Vx**2*qAi + 2*KAB*Q*Vx**2*qBi + 2*KAB*Vw**2*cAi*cBi + 2*KAB*Vw*Vx*cAi*qBi + 2*KAB*Vw*Vx*cBi*qAi + 2*KAB*Vx**2*qAi*qBi + Q**2*Vx**2 - 2*Q*Vw*Vx*cBi - 2*Q*Vx**2*qBi + Vw**2*cBi**2 + 2*Vw*Vx*cBi*qBi + Vx**2*qBi**2))/(2*Vw*(KAB - 1))
    qBf1 = (KAB*Q*Vx - KAB*Vw*cAi - KAB*Vx*qAi - Q*Vx - Vw*cBi - Vx*qBi - sqrt(KAB**2*Q**2*Vx**2 - 2*KAB**2*Q*Vw*Vx*cAi - 2*KAB**2*Q*Vx**2*qAi + KAB**2*Vw**2*cAi**2 + 2*KAB**2*Vw*Vx*cAi*qAi + KAB**2*Vx**2*qAi**2 - 2*KAB*Q**2*Vx**2 + 2*KAB*Q*Vw*Vx*cAi + 2*KAB*Q*Vw*Vx*cBi + 2*KAB*Q*Vx**2*qAi + 2*KAB*Q*Vx**2*qBi + 2*KAB*Vw**2*cAi*cBi + 2*KAB*Vw*Vx*cAi*qBi + 2*KAB*Vw*Vx*cBi*qAi + 2*KAB*Vx**2*qAi*qBi + Q**2*Vx**2 - 2*Q*Vw*Vx*cBi - 2*Q*Vx**2*qBi + Vw**2*cBi**2 + 2*Vw*Vx*cBi*qBi + Vx**2*qBi**2))/(2*Vx*(KAB - 1))
      
    result1 = np.array([cAf1, cBf1, qAf1, qBf1])
    
    ######## solution 2 ##########
    
    cAf2 = (-KAB*Q*Vx + KAB*Vw*cAi + KAB*Vx*qAi + Q*Vx - 2*Vw*cAi - Vw*cBi - 2*Vx*qAi - Vx*qBi + sqrt(KAB**2*Q**2*Vx**2 - 2*KAB**2*Q*Vw*Vx*cAi - 2*KAB**2*Q*Vx**2*qAi + KAB**2*Vw**2*cAi**2 + 2*KAB**2*Vw*Vx*cAi*qAi + KAB**2*Vx**2*qAi**2 - 2*KAB*Q**2*Vx**2 + 2*KAB*Q*Vw*Vx*cAi + 2*KAB*Q*Vw*Vx*cBi + 2*KAB*Q*Vx**2*qAi + 2*KAB*Q*Vx**2*qBi + 2*KAB*Vw**2*cAi*cBi + 2*KAB*Vw*Vx*cAi*qBi + 2*KAB*Vw*Vx*cBi*qAi + 2*KAB*Vx**2*qAi*qBi + Q**2*Vx**2 - 2*Q*Vw*Vx*cBi - 2*Q*Vx**2*qBi + Vw**2*cBi**2 + 2*Vw*Vx*cBi*qBi + Vx**2*qBi**2))/(2*Vw*(KAB - 1))
    qAf2 = (KAB*Q*Vx + KAB*Vw*cAi + KAB*Vx*qAi - Q*Vx + Vw*cBi + Vx*qBi - sqrt(KAB**2*Q**2*Vx**2 - 2*KAB**2*Q*Vw*Vx*cAi - 2*KAB**2*Q*Vx**2*qAi + KAB**2*Vw**2*cAi**2 + 2*KAB**2*Vw*Vx*cAi*qAi + KAB**2*Vx**2*qAi**2 - 2*KAB*Q**2*Vx**2 + 2*KAB*Q*Vw*Vx*cAi + 2*KAB*Q*Vw*Vx*cBi + 2*KAB*Q*Vx**2*qAi + 2*KAB*Q*Vx**2*qBi + 2*KAB*Vw**2*cAi*cBi + 2*KAB*Vw*Vx*cAi*qBi + 2*KAB*Vw*Vx*cBi*qAi + 2*KAB*Vx**2*qAi*qBi + Q**2*Vx**2 - 2*Q*Vw*Vx*cBi - 2*Q*Vx**2*qBi + Vw**2*cBi**2 + 2*Vw*Vx*cBi*qBi + Vx**2*qBi**2))/(2*Vx*(KAB - 1))
    cBf2 = (-KAB*Q*Vx + KAB*Vw*cAi + 2*KAB*Vw*cBi + KAB*Vx*qAi + 2*KAB*Vx*qBi + Q*Vx - Vw*cBi - Vx*qBi - sqrt(KAB**2*Q**2*Vx**2 - 2*KAB**2*Q*Vw*Vx*cAi - 2*KAB**2*Q*Vx**2*qAi + KAB**2*Vw**2*cAi**2 + 2*KAB**2*Vw*Vx*cAi*qAi + KAB**2*Vx**2*qAi**2 - 2*KAB*Q**2*Vx**2 + 2*KAB*Q*Vw*Vx*cAi + 2*KAB*Q*Vw*Vx*cBi + 2*KAB*Q*Vx**2*qAi + 2*KAB*Q*Vx**2*qBi + 2*KAB*Vw**2*cAi*cBi + 2*KAB*Vw*Vx*cAi*qBi + 2*KAB*Vw*Vx*cBi*qAi + 2*KAB*Vx**2*qAi*qBi + Q**2*Vx**2 - 2*Q*Vw*Vx*cBi - 2*Q*Vx**2*qBi + Vw**2*cBi**2 + 2*Vw*Vx*cBi*qBi + Vx**2*qBi**2))/(2*Vw*(KAB - 1))
    qBf2 = (KAB*Q*Vx - KAB*Vw*cAi - KAB*Vx*qAi - Q*Vx - Vw*cBi - Vx*qBi + sqrt(KAB**2*Q**2*Vx**2 - 2*KAB**2*Q*Vw*Vx*cAi - 2*KAB**2*Q*Vx**2*qAi + KAB**2*Vw**2*cAi**2 + 2*KAB**2*Vw*Vx*cAi*qAi + KAB**2*Vx**2*qAi**2 - 2*KAB*Q**2*Vx**2 + 2*KAB*Q*Vw*Vx*cAi + 2*KAB*Q*Vw*Vx*cBi + 2*KAB*Q*Vx**2*qAi + 2*KAB*Q*Vx**2*qBi + 2*KAB*Vw**2*cAi*cBi + 2*KAB*Vw*Vx*cAi*qBi + 2*KAB*Vw*Vx*cBi*qAi + 2*KAB*Vx**2*qAi*qBi + Q**2*Vx**2 - 2*Q*Vw*Vx*cBi - 2*Q*Vx**2*qBi + Vw**2*cBi**2 + 2*Vw*Vx*cBi*qBi + Vx**2*qBi**2))/(2*Vx*(KAB - 1))
        
    result2 = np.array([cAf2, cBf2, qAf2, qBf2])
    
    if (result1<0).any():
        result = result2
    else:
        result = result1
    
    if (result<0).any():
        print('Negative concentations. Something went wrong')
        
    return result


class HSDMIXTestSuite(unittest.TestCase):
    """Test HSDMIX model."""
    
    def setUp(self):
        warnings.simplefilter('ignore', PendingDeprecationWarning)

    def test_trace_film(self):
        """ Film Batch Test """
        
        IEX2 = hsdmixbatch.HSDMIXbatch(DATA_DIR+'/trace_test_hsdmixbatch.xlsx')
        t2, u2 = IEX2.solve()
        
        kL = IEX2.params.kL[0]   # liquid film transfer coefficient (cm/s)
        KAB = IEX2.ions.iloc[1,1] # Selectivity coefficient
        rb = IEX2.params.rb[0] # resin bead radius (cm)
        
        Vw = IEX2.params.VL[0] # volume of water (L)
        Vx = IEX2.params.VR[0] # volume of polymer (L)
        
        Q = IEX2.params.Qm * IEX2.params.RHOP # resin capacity meq/L of wet resin phase
        
        cAi = 0.01 # Contaminant
        cBi = 5  # presat in solution
        qAi = 0 # resin initially free of contaminant
        qBi = Q # resin fully loaded with presat
                
        result = binary_IEX_eq(cAi, cBi, qAi, qBi, Vw, Vx, KAB)   
        cAf, cBf, qAf, qBf = result
        
        tvals = t2

        C_C0 = helfferich_film(tvals, Q, KAB, rb, kL, cAi, cBi, Vw, Vx)
        concs = C_C0 * cAi
        
        assert np.allclose(concs, u2[0, -1, :], rtol=0.01)

if __name__ == '__main__':
    unittest.main()
