# -*- coding: utf-8 -*-
"""
Homogeneous Surface Diffusion Model (HSDM) for ion exchange (IX)

For an overview of ion exchange column modeling, consult:
  Slater, M.J., 2013. Principles of ion exchange technology. Butterworth-Heinemann.
and
  Helfferich, F. G. (1995). Ion exchange. Courier Corporation.

For details on the numerical method of solution (Orthogonal Collocation), consult:
  Crittenden, J. C., Hutzler, N. J., Geyer, D. G., Oravitz, J. L., & Friedman, G. (1986). 
  Transport of organic compounds with saturated groundwater flow: Model development and 
  parameter sensitivity. Water Resources Research, 22(3), 271-284.
  
Assumptions:
  Constant selectivity.
  Plug flow.
  Fickian diffusion.
  Common mass transport parameters for all species.

TODO: Option to calculate film transfer coefficient from correlation.

TODO: Clean, test, double check bicarb/alka output

XXX: Needs a way to specify max_step in solve_ivp to avoid missing influent features
     (NOTE: t_eval doesn't solve this problem. It just interpolates . . .)


@authors: Jonathan Burkhardt, Boris Datsov, Levi Haupert
"""

import timeit

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

from .colloc import build_collocation, advect_operator
from .paramsheets import conv_time, conv_database, conv_conc, conv_params_data

bicarbMW = 61.02

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
    and the order of r is [C, q_r=0, ..., q_r=rb], the ion exchange zones
    will be in the (NION*nz x NION*nz) corners of the Jacobian.
    They have band structure, but we will assume they are dense for simplicity.
    ...
    The Jacobian is dominated by a sparse banded structure from diffusion in the
    solid phase. There are 2*nr - 1 diagonals.
    ...
    Returns Jac_struc: an array of ones and zeros.
    """
    NEQ = (nr+1) * NION * nz
    nzni = NION * nz
    
    Jac_struc = np.zeros((NEQ, NEQ))
    
    # Diffusion zone
    Jac_struc[nzni:, nzni:] += np.eye(NEQ-nzni, NEQ-nzni, k=0)
    for ii in range(1, nr):
        Jac_struc[nzni:, nzni:] += np.eye(NEQ-nzni, NEQ-nzni, k=(ii*nzni))
        Jac_struc[nzni:, nzni:] += np.eye(NEQ-nzni, NEQ-nzni, k=-(ii*nzni))
    
    # Block off corners (ion exchange zones)
    Jac_struc[0:nzni, 0:nzni] = 1.0
    Jac_struc[0:nzni, (nr)*nzni:] = 1.0
    Jac_struc[(nr)*nzni:, 0:nzni] = 1.0
    Jac_struc[(nr)*nzni:, (nr)*nzni:] = 1.0
    
    return Jac_struc


class HSDMIX:
    """ HSDM ion exchange: column process. Plug flow."""
    def __init__(self, inp_file):
        """  """
        self.load_data(inp_file)
        
        
    def load_data(self, inp_file):
        """
        JBB: Unit options. Flow Rate, Diameter. ALKALINITY
        TODO: Deal with case sensitivity: something that capitalizes everything... 
        """
        xls = pd.ExcelFile(inp_file)
        self.params = pd.read_excel(xls, \
                                    sheet_name='params',\
                                    header = [0],\
                                    index_col = [0])
      
        self.params = conv_params_data(self.params)

        self.params = self.params.drop('units', axis=1)['value'] #drops unused column
        
        self.ions = pd.read_excel(xls, \
                                  sheet_name='ions',\
                                  header=[0],\
                                  index_col=[0])
        
        self.valences = self.ions['valence'].values
        
        self.Cin_t = pd.read_excel(xls, \
                                   sheet_name='Cin',\
                                   header=[0], \
                                   index_col = [0],
                                   dtype=np.float64)
        
        
        self.Cin_temp = self.Cin_t.copy(deep=False)
       
        self.time_mult = self.params['time']        
        
        self.Cin_dict = self.ions.to_dict('index')
        
        self.u_Cin2 = {}
        self.val2 = {}
        self.MW2 = {}
        self.u_Cout2 = {}
        
        for c in self.Cin_dict.keys():
            self.u_Cin2[c] = self.Cin_dict[c]['units']
            self.val2[c] = self.Cin_dict[c]['valence']
            self.MW2[c] = self.Cin_dict[c]['mw']
            self.u_Cout2[c] = 'meq'

        self.C_out2, self.u_Cin2, self.u_C_out2 = conv_database(self.Cin_temp, \
                                                               self.u_Cin2, \
                                                               self.u_Cout2, \
                                                               conv_conc, \
                                                               self.MW2, \
                                                               self.val2)
        
        if 'BICARBONATE' not in self.C_out2.columns:
            
            if 'ALKALINITY' in self.C_out2.columns:
#                initialize a column that is the same size as other columns
                self.C_out2['BICARBONATE'] = self.C_out2['ALKALINITY'] * 1.
                self.ions.loc['BICARBONATE'] = self.ions.loc['ALKALINITY']
                self.ions.at['BICARBONATE','mw'] = bicarbMW
                self.ions.at['BICARBONATE','units'] = 'mg'
                
                pH_exp = self.C_out2['PH'] - 10. #convenience
                self.C_out2['BICARBONATE'] = (self.C_out2['ALKALINITY'] - 5. * 10 **pH_exp)/\
                                            (1. + 0.94 * 10**pH_exp)
            else:
                print('Error! No BICARBONATE or ALKALINITY concentration defined')
        
        #clean up un needed columns        
        if 'ALKALINITY' in self.C_out2.columns:
            self.C_out2 = self.C_out2.drop('ALKALINITY', axis=1)
            self.ions = self.ions.drop('ALKALINITY')
        if 'PH' in self.C_out2.columns:
            self.C_out2 = self.C_out2.drop('PH', axis=1)
            
        self.names = self.ions.index.values 
            
    
    def save_results(self, output_file_name, **kwargs):
        '''
        Returns:
            generate and write a *.xlsx file in the parent directory;
            sheet_name = Cout;
            
            *** convert results from solver() to the requested uints; ***
            *** convert results from solver() to the input units if units are not specified; ***
            
        Parameters:
            output_file_name : file name as string;
            period : string;
            units : string; 
            
            *** takes units from the input file if units are not specified; ***
        '''
        
        period = kwargs.get('period', 'hours')
        units = kwargs.get('units', None)
        
        u_Cout2 = {}
        if units == None:
            u_Cout2 = self.u_Cin2
        else:
            for c in self.Cin_dict.keys():
                u_Cout2[c] = units
        
        u_Cin2 = {}
        for c in self.Cin_dict.keys():
            u_Cin2[c] = 'meq'
            
            
        temp_t = pd.Series(self.result.t * self.timeback)
        tmp_u = self.u_result[0,:,-1,:]
        
        
        if period == 'BV':
            bv = temp_t / (self.params['L'] / self.params['v']) 
            idx = pd.Index(bv, name=period)

        else:
            period_factor, u_in, u_out = conv_time('sec',period,'time', period)            
            idx = pd.Index(temp_t * period_factor, name = period)
        
        df_c = pd.DataFrame(tmp_u.T, index = idx, columns = self.names)
        
        tmp_df = df_c.copy(deep=True)

        C_out3, u_Cin2, u_Cout2 = conv_database(tmp_df, u_Cin2, u_Cout2, \
                                    conv_conc, self.MW2, self.val2)       
        
        col_name = C_out3.columns.tolist()
        
        new_col_name = []
        for n in col_name:
            new_col_name.append( n + ' ('+ u_Cout2[n] + '/L)')           

        C_out3.columns = new_col_name
        saved_name = output_file_name
        C_out3.to_excel(output_file_name, sheet_name = 'Cout', float_format='%.8f')        

        return saved_name
    
    
    def solve(self, t_eval=None, const_Cin=False, OCFE=False, quiet=True):
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
        L = np.float64(self.params['L']) # Column Length (cm)
        v = np.float64(self.params['v']) # linear flow velocity (cm/s)
        nz = int(self.params['nz']) # number of axial collocation points
        nr = int(self.params['nr']) # number of radial collocation points
        Ds = np.float64(self.params['Ds']) # overall resin phase diffusion coefficient (cm**2/s)
        kL =  np.float64(self.params['kL'])   # overall film transfer coefficient (cm/s)
        rb = np.float64(self.params['rb']) # resin radius (cm)
        Q = np.float64(self.params['Q']) # Resin capacity by volume (meq/L)

        
        ##########################
        ### DERIVED PARAMETERS ###
        ##########################
        
        CT = Cin.sum() # XXX: Need to modify for time varying feed concentrations
        DGT = (1 - EBED) / EBED * Q / CT  # overall distribution coefficient :XXXX:
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
        u[RESIN:, PRESAT, :] = Q  # resin initially loaded with PRESAT

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
            Cc_eq = (-bb + np.sqrt(bb**2 - 4 * aa * cc))/(2*aa)

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
            u[RESIN:, dv_idx, :] = Q/1e3 # so the quadratic formula doesn't explode
            u[RESIN:, PRESAT, :] = Q - (ndv * Q / 1e3)  
            
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
            q = u[RESIN:, :, :]
            q_s = u[SURF, :, :]
            
            u[LIQUID, :, 0] = get_Cin(T)    # XXX: SLOW
            
            # Calculate total concentrations along column
            CT = u[LIQUID, :, :].sum(axis=0)
            
            # update Ceq
            Ceq = calc_Ceq(q_s, CT)
            
            # Calculate flux terms
            J = - kL * (C - Ceq) # mass flux 
            Jas = J * 3/rb  # mass flux * specific surface area of bead

            # Initialize arrays for derivatives
            du_dT = np.zeros([nr+1, NION, nz])
            dq_dT = np.zeros([nr, NION, nz])
            dC_dT = np.zeros([NION, nz]) 
            dq_dT_w = np.zeros([NION, nz])
        
            # Liquid phase
            C_swap = np.swapaxes(C, 0, 1)
            Az_C = np.swapaxes(np.matmul(Az, C_swap), 0, 1)
            dC_dT = (-v/L*(Az_C) + (1-EBED)*Jas)/EBED * t_half
        
            # diffusion in bead  
            q_swap = np.swapaxes(q, 0, 1)
            Br_q = np.swapaxes(np.matmul(Br, q_swap), 0, 1)
            dq_dT =  Ds * t_half / rb**2 * Br_q   

            # intermediate term for dq_dT at bead surface        
            dq_dT_swap = np.swapaxes(dq_dT[:SURF, :, :],0,1)
            dq_dT_w = np.matmul(Wr[:-1], dq_dT_swap)
            
            # Fill out du_dT            
            du_dT[LIQUID, :, :] = dC_dT
            du_dT[LIQUID, :, 0] = 0 # Inlet
            du_dT[RESIN:, :, :] = dq_dT
            du_dT[SURF, :, :] = (-t_half / rb * J - dq_dT_w)/Wr[-1]       
     
            du_dT = du_dT.reshape(NEQ) # solve_ivp requires this to be a vector
           
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
            print('HSDM solve_time (s): ' + str(solve_time))

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
            
        if not np.allclose(u[SURF, :, :, :].sum(axis=0), Q, rtol=0.01):
            print('WARNING: Sum of resin concentrations is not Q!')
            
        if np.any(u[:, :, :, -1] < 0):
            print('WARNING: Negative concentrations detected!')      
        
        return (t, u)


