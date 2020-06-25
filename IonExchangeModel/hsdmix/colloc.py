# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 12:18:19 2018

Generate Collocation operators.

Reference: Finlayson's "numerical methods for chemical engineers"

@author: LHaupert
"""

import numpy as np
from scipy.special import roots_sh_jacobi
from scipy.special import roots_sh_legendre
from numpy.linalg import inv

def build_collocation(nr, nz):
    """
    nr = number of radial points
    nz = number of axial points
    
    Returns
    Az = axial first derivative operator
    Br = 1D radial Laplacian operator for a sphere
    Wr = vector of radial quadrature weights
    """

    #############################################
    # Make Axial first derivative operator Az ###
    #############################################
    n_pts = nz # total number of collocation points
    n_pts_m = n_pts - 2 # number of interior collocation points
    
    rootsz = np.zeros([n_pts])
    rootsz[0] = 0    # left boundary root
    rootsz[-1] = 1   # right bounary root
    rootsz[1:-1] = roots_sh_legendre(n_pts_m)[0]
    
    Qmat = np.zeros([n_pts, n_pts])
    Cmat = np.zeros([n_pts, n_pts])
    
    for i in range(n_pts):
        for j in range(n_pts):
            Qmat[i, j] = rootsz[i]**j
    
        for j in range(1, n_pts):
            Cmat[i, j] = (j)*rootsz[i]**[j-1]
    
    Qinv = inv(Qmat)
    Az = np.dot(Cmat, Qinv)
    
    ###############################################################
    # Make sphere Laplacian Br, and quatrature weights Wr         #
    ###############################################################
    
    n_pts = nr # total number of collocation points
    n_pts_m = n_pts - 1 # number of interior collocation points
    # NOTE symmetric around x = 0
    
    rootsr = np.zeros([n_pts])
    rootsr[-1] = 1.   # right bounary root
    rootsr[:-1] = np.sqrt(roots_sh_jacobi(n_pts_m, 2.5, 1.5)[0])
    
    Qmat = np.zeros([n_pts, n_pts])
    Dmat = np.zeros([n_pts, n_pts])
    
    for i in range(n_pts):
        for j in range(n_pts):
            Qmat[i, j] = rootsr[i]**(2*j)
    
        for j in range(1, n_pts):
            Dmat[i, j] = (2*j)*(2*j+1)*rootsr[i]**[2*j-2]  # NOTE: spherical symmetry
    
    Qinv = inv(Qmat)
    Br = np.dot(Dmat, Qinv)
    
    # Calculate quadrature weights
    # (See Finlayson's "numerical methods for chemical engineers")
    ii = np.arange(n_pts) + 1
    f = 1./(2*(ii) -2 + 3)
    Wr = np.dot(f, Qinv)    # quadrature weighting vector
    
    return rootsz, Az, rootsr, Br, Wr


def advect_operator(NE, NP):
    """
    NE: Number of elements of uniform with
    NP: Number of collocation points per element.
    
    Returns:
        Advection operator
        location of axial collocation points
    """
    
    NZ = NE * NP - (NE-1) # Number of axial points
    F = NE # Element width adjustment to derivatives
    
    roots, A, _, _, _ = build_collocation(3, NP)
    
    
    ### Calculation location of overall gridpoints.
    Xvals = np.zeros(NE * NP)
    for k in range(NE):
        Xvals[k*NP:(k+1)*NP] = roots/F + k/F
        
    to_delete = [NP*(k+1) for k in range(NE-1)]
    Xvals = np.delete(Xvals, to_delete)  
    
    Q = A[-1, -1] - A[0, 0] # For convenience later

    ## Adjusted element advection operator
    Z = np.zeros([NP, NP])
    Z[:, :] = F*A

    ###  Construct overall advection operator
    Adv_Op = np.zeros([NZ, NZ])

    # Fill in blocks for element interiors
    for k in range(NE):
        ii = k*NP-k
        iii = (k+1)*NP-k
        Adv_Op[ii:iii, ii:iii] = Z[:, :]
        
   
    # Fill in continuation points where elements meet
    for k in range(NE-1):
        idx = k*NP-k
        ii = (k+1)*NP-(k+1)
        iii= (k+2)*NP-(k+1)
        
        Adv_Op[ii, :] = 0  # Zero out continuation
    
        CC1 = Z[-1,:-1] - Z[-1, -1] / Q * A[-1, :-1]
        CC2 = Z[-1, -1] / Q * A[0, 1:]

        Adv_Op[ii, idx:ii] = CC1
        Adv_Op[ii, ii+1:iii] = CC2   # XXX: clearer indexing?

    # Inlet Boundary
    Adv_Op[0, :] = 0   # Constant

        
    return Adv_Op, Xvals    

