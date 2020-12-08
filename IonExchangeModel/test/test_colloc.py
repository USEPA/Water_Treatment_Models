# -*- coding: utf-8 -*-

from context import ixpy

import unittest

import numpy as np


class CollocTestSuite(unittest.TestCase):
    """Test Collocation construction."""


    def test_Br(self):
        _, _, _, Br, _ = ixpy.colloc.build_collocation(3, 4)
        Br_test = np.array([[-15.66996162,  20.03487835,  -4.36491673],
                            [  9.96512165, -44.33003838,  34.36491673],
                            [ 26.93285519, -86.93285519,  60.        ]])
        assert np.allclose(Br_test, Br)
        
        
    def test_Wr(self):
        _, _, _, _, Wr = ixpy.colloc.build_collocation(3, 4)
        Wr_test = np.array([ 0.09490592,  0.19080837,  0.04761905])   
        assert np.allclose(Wr_test, Wr)


    def test_Az(self):
        _, Az, _, _, _ = ixpy.colloc.build_collocation(3, 4)
        Az_test = np.array([[-7.,          8.19615242, -2.19615242,  1.,        ],
                            [-2.73205081,  1.73205081,  1.73205081, -0.73205081,],
                            [ 0.73205081, -1.73205081, -1.73205081,  2.73205081,],
                            [-1.,          2.19615242, -8.19615242,  7.,        ]]) 
        assert np.allclose(Az_test, Az)
        
        
    def test_OCFE_roots(self):
        _, roots_OCFE = ixpy.colloc.advect_operator(2, 3)
        roots_test = np.array([0.  , 0.25, 0.5 , 0.75, 1.  ])
        assert np.allclose(roots_test, roots_OCFE)
        
    
    def test_OCFE_Az(self):
        Az_OCFE, _ = ixpy.colloc.advect_operator(2, 3)
        Az_OCFE_test = np.array([[ 0.,  0.,  0.,  0.,  0.],
                                 [-2.,  0.,  2.,  0.,  0.],
                                 [ 1., -4.,  0.,  4., -1.],
                                 [ 0.,  0., -2.,  0.,  2.],
                                 [ 0.,  0.,  2., -8.,  6.]])
        assert np.allclose(Az_OCFE_test, Az_OCFE)


if __name__ == '__main__':
    unittest.main()
