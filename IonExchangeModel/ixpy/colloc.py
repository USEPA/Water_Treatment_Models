# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 12:18:19 2018

Generate Collocation operators.

Reference:
Michelsen, M. L., & Villadsen, J. (1972). A convenient computational procedure 
for collocation constants. The Chemical Engineering Journal, 4(1), 64-68.

Further references:
Finlayson, B. A. (1980). Nonlinear analysis in chemical engineering.
    
Villadsen, J., & Michelsen, M. L. (1978). Solution of differential equation 
models by polynomial approximation(Book). Englewood Cliffs, N. J., 
Prentice-Hall, Inc., 1978. 460 p.

Villadsen, J. V., & Stewart, W. E. (1967). Solution of boundary-value problems 
by orthogonal collocation. Chemical Engineering Science, 22(11), 1483-1501.

@author: LHaupert
"""

import numpy as np
from scipy.special import roots_sh_jacobi
from scipy.special import roots_sh_legendre


def recur(pj, H):
    """
    Eq. 17
    H = (u - u_j)
    pj = p1_j, p2_j, p3_j
    
    returns p_j+1 (p_next)
    """
    
    p1_next = H * pj[0]
    p2_next = H * pj[1] + 2 * pj[0]
    p3_next = H * pj[2] + 3 * pj[1]
    
    return np.array([p1_next, p2_next, p3_next])


def calc_px_column(root_i, roots):
    """
    Given a root of interest (root_i), and the list of roots, 
    return polynomial derivatives p1, p2, p3 at root_i using eq. 17
    """
    
    idxs = (roots - root_i).nonzero()  # exclude root_i
    
    # initial derivative values for recurrence
    p1 = 1.0
    p2 = 0.0
    p3 = 0.0
    
    pn = np.array([p1, p2, p3])
    
    for jj in idxs[0]:  # do recurrence, root by root, skipping root_i

        H = root_i - roots[jj]
        pn = recur(pn, H)
   
    return pn  # return column for px_mat


def recur_colloc(n_pts):
    """
    For non-symmetric problems produce:
        roots: Root locations
        Amat: First derivative operator
        Bmat: Second derivative 
        
    TODO: Refactor with mv_colloc_symm to avoid code duplication
    """

    n_pts_m = n_pts - 2 # number of interior collocation points
    roots = np.zeros([n_pts])
    roots[0] = 0    # left boundary root
    roots[-1] = 1   # right bounary root
    roots[1:-1] = roots_sh_legendre(n_pts_m)[0] # interior roots
    
    # define 3 by n_pts matrix of polynomial derivatives at collocation points
    # row index (k) is the k + 1 derivative
    px_mat = np.zeros((3, n_pts))
    
    for ii in range(n_pts):
        px_mat[:, ii] = calc_px_column(roots[ii], roots)
    
    # Calculate derivative operators (A and B)
    Amat = np.zeros((n_pts, n_pts))
    Bmat = np.zeros((n_pts, n_pts))
    
    for ii in range(n_pts):
        for jj in range(n_pts):
            if ii == jj:  # Eq. 12
                Amat[ii, jj] = (1/2) * px_mat[1, ii] / px_mat[0, ii]
                Bmat[ii, jj] = (1/3) * px_mat[2, ii] / px_mat[0, ii]
            else: # Eqs. 13, 14
                G = roots[ii] - roots[jj]
                Amat[ii, jj] = (1/G) * px_mat[0, ii] / px_mat[0, jj]
                Bmat[ii, jj] = (1/G) * (px_mat[1, ii] / px_mat[0, jj] - 2 * Amat[ii, jj])
                
    return(roots, Amat, Bmat)


def recur_colloc_symm(n_pts, a):
    """
    For symmetric collocation problems, produce:
        rootsx = polynomial root locations
        Ax = 1st derivative operator
        Lx = 1-D Laplacian
        W = quadrature weights
    Reference: Michelsen and Villadsen 1971 
    sheet: a=1
    cylinder: a=2
    sphere: a=3
    """

    n_pts_m = n_pts - 1 # number of interior collocation points
    rootsu = np.zeros([n_pts])

    rootsu[-1] = 1   # right bounary 
    # Define optimal Jacobi polynomial for interior points
    p = 1.5 + (a-1)/2
    q = p - 1
    rootsu[:-1] = roots_sh_jacobi(n_pts_m, p, q)[0]  # roots in u domain
    rootsx = np.sqrt(rootsu) # roots in x domain
    
    
    # define 3 by n_pts matrix of polynomial derivatives at collocation points
    # row index (k) is the k + 1 derivative
    px_mat = np.zeros((3, n_pts))
    
    for ii in range(n_pts):
        px_mat[:, ii] = calc_px_column(rootsu[ii], rootsu)
    
    # Calculate derivative operators (A and B) in u domain
    Amat = np.zeros((n_pts, n_pts))
    Bmat = np.zeros((n_pts, n_pts))
    
    for ii in range(n_pts):
        for jj in range(n_pts):
            if ii == jj:  # Eq. 12
                Amat[ii, jj] = (1/2) * px_mat[1, ii] / px_mat[0, ii]
                Bmat[ii, jj] = (1/3) * px_mat[2, ii] / px_mat[0, ii]
            else: # Eqs. 13, 14
                G = rootsu[ii] - rootsu[jj]
                Amat[ii, jj] = (1/G) * px_mat[0, ii] / px_mat[0, jj]
                Bmat[ii, jj] = (1/G) * (px_mat[1, ii] / px_mat[0, jj] - 2 * Amat[ii, jj])
                
    Ax = np.zeros((n_pts, n_pts)) # First derivative operator in x domain
    Bx = np.zeros((n_pts, n_pts)) # 1-D Laplacian operator in x domain
    for ii in range(n_pts):
        for jj in range(n_pts):
            # Following Eq. 33
            Ax[ii, jj] = 2 * rootsx[ii] * Amat[ii, jj]
            Bx[ii, jj] = 2 * a * Amat[ii, jj] + rootsu[ii] * 4 * Bmat[ii, jj]
            
    W = np.zeros(n_pts)
    for ii in range(n_pts):
        # Eq. 23 - 25 but omitting root at zero for nodal polynomial.
        wp_sum = (1 / (rootsu[:]*px_mat[0, :]**2)).sum()
        W[ii] = (1 / (rootsu[ii] * px_mat[0, ii]**2)) / wp_sum / (a)
    
    return(rootsx, Ax, Bx, W)


def build_collocation(nr, nz):
    """
    nr = number of radial points
    nz = number of axial points
    
    Returns
    rootsz = axial roots
    Az = axial first derivative operator
    rootsr = radial roots
    Br = 1D radial Laplacian operator for a sphere
    Wr = vector of radial quadrature weights
    """

    rootsr, _, Br, Wr = recur_colloc_symm(nr, 3)
    rootsz, Az, _ = recur_colloc(nz)
    
    return rootsz, Az, rootsr, Br, Wr


def advect_operator(NE, NP):
    """
    For axial OCFE
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

