# -*- coding: utf-8 -*-
"""
Pore-Surface Diffusion Model (PSDM) for ion exchange (IX)

Extends HSDMIX. Consult hsdmix.py for details
     
TODO: Use further inheritance or composition (or something) to reduce code 
      duplication with hsdmix


@authors: Jonathan Burkhardt, Boris Datsov, Levi Haupert
"""

import timeit
import pandas as pd
import numpy as np
import copy
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

from .colloc import build_collocation, advect_operator

from .hsdmix import HSDMIX
from .converter import conv_iex_u

def approx_Jac_struc(nr, NION, nz):
    """
    Find approximate Jacobian structure to speed up BDF or Radau calculations
    in solve_ivp
    ...
    nr: Number of radial collocation points
    NION: Number of species
    nz: Number of axial collocation points
    ...
    Because the order of the concentration array is (nr, NION, nz),
    and the order of r is [C, q_r=0, ..., q_r=rb].
    ...
    Returns Jac_struc: an array of ones and zeros.
    """
    NEQ = (nr+1) * NION * nz
    nzni = NION * nz
    
    Jac_struc = np.zeros((NEQ, NEQ))
    
    # Diffusion/exchange zone
    Jac_struc[nzni:, nzni:] += np.eye(NEQ-nzni, NEQ-nzni, k=0)
    
    # No axial interactions in beads
    for ii in range(NION*(1+nr)-1):
        ii = ii+1
      
        # just do the whole block and zero out later. It's easier.
        Jac_struc[:, :] += np.eye(NEQ, NEQ, k=nz*ii)
        Jac_struc[:, :] += np.eye(NEQ, NEQ, k=-nz*ii)
    
    # zero out liquid phase
    Jac_struc[0:nzni,:] = 0.0        
    
    # Block off corners (ion exchange zones)
    Jac_struc[0:nzni, 0:nzni] = 1.0 # advection zones
    Jac_struc[0:nzni, (nr)*nzni:] = 1.0 # bead surface
    
    return Jac_struc


class PSDMIX(HSDMIX):
    """ PSDM ion exchange: column process. Plug flow."""
    def __init__(self, inp_file):
        """  """
        super().__init__(inp_file)  # Use HSDMIX input processing machinery
        # Note: also inherits save_results method from HSDMIX
    
    
    def solve(self, t_eval=None, const_Cin=False, OCFE=False, quiet=True, u_init=None):
        """ Returns (t, u)
        t = time values
        u = array of shape (phases, ions, axial, time) 
        Liquid phase = 0 
        Resin phase = 1:
        t_eval = requested reporting timesteps
        const_Cin: if True, constant influent concentrations are assumed
        OCFE: if True, orthogonal collocation on finite elements will be used in z
              note that nz must be an odd number.
        """
        
        Cin = [] # initial inlet concentrations
        for name in self.names:
            Cin.append(self.Cin_t[name].values[0])
        Cin = np.array(Cin)
        self.Cin = Cin
        
        
        ### Alias some parameters for convenience ###

        EBED = np.float64(self.params['EBED']) # bed porosity
        EPOR = np.float64(self.params['EPOR']) # bead porosity
        L = np.float64(self.params['L']) # Column Length (cm)
        v = np.float64(self.params['v']) # linear flow velocity (cm/s)
        nz = int(self.params['nz']) # number of axial collocation points
        nr = int(self.params['nr']) # number of radial collocation points
        Ds = np.float64(self.params['Ds']) # overall resin phase diffusion coefficient (cm**2/s)
        Dp = np.float64(self.params['Dp']) # pore diffusion coefficient (cm**2/s)
        kL =  np.float64(self.params['kL'])   # overall film transfer coefficient (cm/s)
        rb = np.float64(self.params['rb']) # resin radius (cm)
        Y = np.float64(self.params['Q']) # Resin capacity by bead volume (meq/L)

        
        ##########################
        ### DERIVED PARAMETERS ###
        ##########################
        
        
        CT = Cin.sum() # XXX: Need to modify for time varying feed concentrations
        DGT = (1 - EBED) / EBED * Y / CT  # overall distribution coefficient :XXXX:
        tau = L * EBED / v # tracer residence time
        t_half =  tau * (1 + DGT) # time for outlet conc to be half presat
        self.timeback = t_half
        # NOTE: T = t / t_half
        self.toBV = EBED / tau # May be useful for analyzing data
        
        NION = len(Cin) 
        NEQ = (nr + 1) * (NION) * nz # number of equations 
        
        Jac_sp = approx_Jac_struc(nr, NION, nz) # approximate structure of jacobian
        
        # Enumeration to keep track of indexes. 
        LIQUID = 0
        RESIN = 1
        PRESAT = 0
        SURF = -1
        
        # divalent indexes
        valences = self.valences # XXX: Consider cleaning this up
        val2 = np.tile(np.array([valences]).T, (1, nz))   # tiled for convenience
        dv_idx = valences == 2
        ndv = len(valences[dv_idx])# number of divalent species
        
        # non chloride monovalent indexes
        mv_idx = valences == 1
        mv_idx[0] = False   # presaturant (reference ion) 
        
        # Equilibrium coefficients
        Kxc = self.ions['Kxc'].values  # separation factor vs presat, 
        # NOTE: the FIRST entry is the presaturant
        Kxc2 = np.tile(np.array([Kxc]).T, (1, nz)) # Tiled for convenience

        # MONOVALENT ONLY! Construct separation factor matrix.
        alpha_in = self.ions['Kxc'].values  # separation factor vs presat, 
        # NOTE: the FIRST entry in alpha_in is the presaturant
        alpha = np.zeros([NION, NION])
        for i in range(NION):
            for j in range(NION):
                alpha[i, j] = alpha_in[i]/alpha_in[j]  
        
        

        ##############################
        ### ORTHOGONAL COLLOCATION ###
        ##############################

        rootsz, Az, rootsr, Br, Wr = build_collocation(nr, nz)
        if OCFE:
            if not nz % 2: # nz is even
                print('WARNING: nz must be an odd number in OCFE mode.')
            NE = int((nz-1)/2)
            Az, rootsz = advect_operator(NE, 3) 
        # Az is the spatial first derivative operator
        # Br is the radial, symmetric second derivative operator for a sphere
        # Wr is a vector of quadrature weights
        # rootsz and rootsr are the coordinates of the collocation points

        ### Initialize concentrations ###
        u = np.zeros([(nr + 1), NION, nz]) # equivalent concentration tensor
        for i in range(NION):
            u[LIQUID, i, 0] = Cin[i]
        
        u[LIQUID, PRESAT, 1:] = Cin.sum() # column initially full of presat solution?
        u[RESIN:, PRESAT, :] = Y  # resin initially loaded with PRESAT
        if np.any(u_init): # start with a different distribution in the resin.
            u[RESIN:, :, :] = u_init[RESIN:, :, :]
             
        u0 = u.reshape(NEQ)  # initial concentration vector for solve_ivp call

        
        #######################
        ### Time Series Cin ###
        #######################
        
        
        interp_list = []
        
        # XXX: correct way to deal w/ time units?
        to_T = self.time_mult / t_half ### handles time? but then messes up lines below JBB
        
        for name in self.names:
            #make time interpolating function
            # XXX: What if names don't match?!
            finterp = interp1d(self.Cin_t.index * to_T, self.Cin_t[name], #### self.time_mult messes this up???? JBB
                                          kind='linear')
            interp_list.append(finterp)
            
        
        def get_Cin(T, const=const_Cin):
            # XXX: SLOW! Need to vectorize or otherwise speed up somehow.
            # Doing this almost quadruples solve time.
            # But I don't think this can be done w/o Cython, if possible at all.
            if const:
                return self.Cin
            else:
                Cin_list = [f(T) for f in interp_list]
                return np.array(Cin_list)
        
        self.get_Cin = get_Cin  # This will allow us to access influent interpolation
        # XXX: Should probably move some of these function definitions
        # out of solve

        #########################
        ### Local Equilibrium ###
        #########################
        
        def calc_Ceq_dv(q_s, CT):
            # update Ceq with proper accounting for divalent ions
            # XXX: Slower than calculating using separation factor.
            # Appears to cause additional instability in advection part
            # There are actually several ways this could go wrong.
            # XXX: What if q_divalent drops to zero during the run?
            Ceq = np.zeros(np.array([NION, nz]))
            
            aa = (q_s[dv_idx, :] / Kxc2[dv_idx, :]).sum(axis=0) / q_s[PRESAT, :]**2
            bb = (q_s[mv_idx, :] / Kxc2[mv_idx, :]).sum(axis=0) / q_s[PRESAT, :] + 1
            cc = -CT
   
            Cc_eq = np.zeros([nz])
            Cc_eq = 2 * cc / (-bb - np.sqrt(bb**2 - 4 * aa * cc))  

            Cc_tile = np.tile(Cc_eq, (NION, 1))
            
            Ceq[0:, :] = q_s[0:, :]/Kxc2[0:, :]*(Cc_tile[0:, :]/q_s[0, :])**val2[0:, :]
                
            Ceq[PRESAT, :] = Cc_eq
            
            return Ceq
        
   
        def calc_Ceq_mv(q, CT):
            denoms = np.dot(alpha, q) # a vector of denominators
            Ceq = q * CT / denoms     # elementwise division 
            return Ceq
        
        if np.any(valences == 2):
            calc_Ceq = calc_Ceq_dv
            
        else:
            print('No divalent ions in input. Interpreting Kxc as separation factor.')
            calc_Ceq = calc_Ceq_mv
        
        
        ###########################
        ### Solving Derivatives ###
        ###########################
        
        def diffun(T, u):
            
            """
            Calculate time derivatives at grid points
            """
            u = u.reshape([(nr + 1), NION, nz])
            
            C = u[LIQUID, :, :]
            Y = u[RESIN:, :, :]
            Y_s = u[SURF, :, :]
            Cpore = np.zeros(Y.shape)
            
            q = Y / (1 - EPOR) # Approximation
            q_s = Y_s / (1 - EPOR) # Approximation
            
            u[LIQUID, :, 0] = get_Cin(T)    # XXX: SLOW
            
            # Calculate total concentrations along column
            CT = u[LIQUID, :, :].sum(axis=0)

            
            # update Ceq
            Ceq = calc_Ceq(q_s, CT)
            
            for ii in range(nr):   # XXX: Slow, but something
                Cpore[ii, :, :] = calc_Ceq(q[ii, :, :], CT)
            
            # Calculate flux terms
            J = np.zeros((NION, nz))
            for iii in range(NION):
                J[iii, :] = - self.ions['kL'][iii] * (C[iii,:] - Ceq[iii,:]) # mass flux 
                
            # explicitly doing implicit chloride
            J[0, :] = -J[1:, :].sum(axis=0)
            Jas = J * 3/rb  # mass flux * specific surface area of bead

            # Initialize arrays for derivatives
            du_dT = np.zeros([nr+1, NION, nz])
            dY_dT = np.zeros([nr, NION, nz])
            dC_dT = np.zeros([NION, nz]) 
            dY_dT_w = np.zeros([NION, nz])
        
            # Liquid phase
            C_swap = np.swapaxes(C, 0, 1)
            Az_C = np.swapaxes(np.matmul(Az, C_swap), 0, 1)
            dC_dT = (-v/L*(Az_C) + (1-EBED)*Jas)/EBED * t_half
            
           
            # diffusion in bead  
            Y_swap = np.swapaxes(Y, 0, 1)
            Br_Y = np.swapaxes(np.matmul(Br, Y_swap), 0, 1)

            
            # pore diffusion
            Cpore_swap = np.swapaxes(Cpore, 0, 1)
            Br_Cpore =  np.swapaxes(np.matmul(Br, Cpore_swap), 0, 1)

            
            for iii in range(NION):
                Dp_iii = self.ions['Dp'][iii]
                Ds_iii = self.ions['Ds'][iii]
                dY_dT[:, iii, :] =  t_half * (EPOR * (Dp_iii - Ds_iii) * Br_Cpore[:, iii, :] + Ds_iii * Br_Y[:, iii, :]) / rb**2

            # # explicitly doing implicit chloride
            # dq_dT[:, 0, :] = -dq_dT[:, 1:, :].sum(axis=1) # XXX: Why doesn't work?
            # print(dq_dT[-2, :, -1].sum()) # Why isn't that zero?

            # intermediate term for dq_dT at bead surface        
            dY_dT_swap = np.swapaxes(dY_dT[:SURF, :, :], 0, 1)
            dY_dT_w = np.matmul(Wr[:-1], dY_dT_swap)

            
            # Fill out du_dT            
            du_dT[LIQUID, :, :] = dC_dT.reshape([NION, nz])
            du_dT[LIQUID, :, 0] = 0 # Inlet
            du_dT[RESIN:, :, :] = dY_dT
            du_dT[SURF, :, :] = (-t_half / rb * J - dY_dT_w)/Wr[-1]       
     
            du_dT = du_dT.reshape(NEQ) # solve_ivp requires this to be a vector
            
            if np.any(np.isnan(du_dT)):
                # XXX: Doesn't solve_ivp already check for this?
                raise ValueError('###### WARNING: At least one derivative is NAN!')
           
            return du_dT  
        
        
        #####################
        ### Solve Problem ###
        #####################
        
        T_final = self.Cin_t.index.values[-1] * to_T ##### changed JBB
        Tvals = np.array([0, T_final]) # time points to solve

        T_eval = None
        if np.any(t_eval):  # specific times requested
            T_eval = t_eval * to_T
        
        start_time = timeit.default_timer()
        self.result = solve_ivp(diffun, Tvals, u0, method='BDF',
                                t_eval=T_eval, jac_sparsity=Jac_sp) 
        solve_time = timeit.default_timer() - start_time
        if not quiet:
            print('PSDM solve_time (s): ' + str(solve_time))

        t = self.result.t * t_half # convert back to actual time
        NT = len(self.result.t)
        u = self.result.y.reshape([(nr + 1), NION, nz, NT])
        self.u_result = u
        
        ########################
        ### Check for Errors ###
        ########################
           
        CT_out = u[0, :, -1, :].sum(axis=0)
        Cin_check = get_Cin(self.result.t)
        CT_in = np.array(Cin_check).sum(axis=0)

        if not np.allclose(CT_out, CT_in, rtol=0.02):  # XXX: Is this tight enough?
            print('WARNING: Total outlet does not follow total inlet!')
            
        if not np.allclose(u[SURF, :, :, :].sum(axis=0), Y, rtol=0.01):
            print('WARNING: Sum of resin concentrations is not Y!')
            
        if np.any(u[:, :, :, -1] < 0):
            print('WARNING: Negative concentrations detected!')      
        
        return (t, u)

    def model_uncertainty_PSDM(self, resin_capacity=10, Ds='None', kL='None', flrt='None', c0='None', Kxc='None', L='None', Dp='None', ebed='None', epor='None'):
        obj_copy = copy.deepcopy(self)

        duration = self.Cin_t.index[-1]
        
        ## adjust tvals for evaluations based on the time multiplier. Reset things to hours.
        if self.time_mult == 3600:
            tvals = np.arange(0, duration+1, 12) ## currently set at every hour
        else: ## assumes days
            self.Cin_t.index *= 24
            tvals = np.arange(0, duration*24 + 1, 12)
            self.time_mult = 3600.
            self.params.loc['time'] = 3600.
        
        ## set up output dataframe that stores results
        m_idx = [(i, j) for i in ['base', 'upper', 'lower'] for j in self.ions.index ]
        midx = pd.MultiIndex.from_tuples(m_idx)
        output_df = pd.DataFrame(columns=midx, index=np.round(tvals/24, 2))
        
        ## run base case
        t, u = self.solve(t_eval=tvals, const_Cin=False, OCFE=False)
        converted_df = conv_iex_u(u, t, self.ions)

        output_df['base'] = converted_df
        
        upper_df = converted_df.values * 1.0
        lower_df = converted_df.values * 1.0

        inputs = {'capacity': resin_capacity, 'Kxc': Kxc, 'ds': Ds, 'kL': kL, 'dp': Dp,
                  'flrt': flrt, 'c0': c0, 'L': L, 'ebed': ebed, 'epor': epor} ## TODO: add porosity/EBED or similar for PSDM? Dp for PSDM

        test_uncertainty = []
        for key, value in inputs.items():
            if value != 'None':
                test_uncertainty.append({key: 1 + value/100})
                test_uncertainty.append({key: 1 - value/100})

        ## TODO: add a simulation where all uncertainties are calculated together???

        for test in test_uncertainty:
            self = copy.deepcopy(obj_copy)  ## reset column

            for key, value in test.items(): ## should be a dictionary
                if key == 'capacity': #Resin Capacity
                    # print(self.params)
                    self.params.loc['Qf'] *= value
                    self.params.loc['Q'] *= value ## not sure if we need to update just Q or just Qf, or both... TODO: Check with Levi
                if key == 'ebed': #Bed Porosity
                    self.params.loc['EBED'] *= value    
                if key == 'epor': ## BEAD Porosity
                    self.params.loc['EBED'] *= value
                if key == 'ds': #Surface Diffusion Coefficient
                    self.params.loc['Ds'] *= value
                if key == 'dp': #Pore Diffusion Coefficient
                    self.params.loc['Dp'] *= value
                if key == 'kL': #Film Transfer Coeffient
                    self.params.loc['kL'] *= value
                if key == 'flrt': # Flow rate
                    ## since velocity is calculated, it can be used directly
                    ## TODO: Check if flrt and or diam need to be recalculated. Don't think these are used in the actual function.
                    self.params.loc['v'] *= value
                if key == 'L': # Bed length
                    self.params.loc['L'] *= value
                if key == 'c0': # Influent Concentrations
                    self.Cin_t *= value 
                if key == 'Kxc': # Selectivity coefficients
                    self.ions['Kxc'] *= value
                    self.ions.loc['CHLORIDE', 'Kxc'] = 1.0 ## reset Chloride KxA to 1, as reference ion


            ### run new simulation with updated parameters
            t, u = self.solve(t_eval=tvals, const_Cin=False, OCFE=False)
            intermediate_df = conv_iex_u(u, t, self.ions)

            ## find values that increase or decrease uncertainty ranges
            upper_df[np.where(upper_df < intermediate_df.values)] = intermediate_df.values[np.where(upper_df < intermediate_df.values)]
            lower_df[np.where(lower_df > intermediate_df.values)] = intermediate_df.values[np.where(lower_df > intermediate_df.values)]

        if len(test_uncertainty) > 0:
            output_df['upper'] = upper_df.copy()
            output_df['lower'] = lower_df.copy()




        self = copy.deepcopy(obj_copy) ## reset class object to original state


        return output_df
