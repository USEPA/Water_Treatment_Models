# -*- coding: utf-8 -*-
"""
Created on Tue May 14 08:24:46 2019

HSDM batch model

@author: LHaupert
"""

import numpy as np
import pandas as pd
from scipy.linalg import block_diag
from scipy.integrate import solve_ivp

from .paramsheets import conv_params_data, conv_database, conv_conc, conv_time
from .colloc import build_collocation


class HSDMIXbatch:
    """ HSDM ion exchange: column process. Plug flow."""
    def __init__(self, inp_file):
        """  """
        self.load_data(inp_file)
        
    def load_data(self, inp_file):
        """ Load parameters, compound data, and concetration series"""
        xls = pd.ExcelFile(inp_file)
        self.params = pd.read_excel(xls, 'params')
        self.params.set_index('name', inplace=True)
        self.paramunits = self.params['units'].copy()
        self.params = conv_params_data(self.params)
        self.params.drop(columns=['units'], inplace=True)
        self.params = self.params.transpose()
        self.ions = pd.read_excel(xls, 'ions')
      
        self.names = self.ions['name'].values    
        self.valences = self.ions['valence'].values
        self.ions.set_index('name', inplace=True)

        self.ions['C0'], _, _ = conv_database(self.ions['C0'], 
                                        self.ions['units'].to_dict(),
                                        dict(zip(self.ions.index.tolist(), ['meq']*len(self.ions['C0']))),
                                        conv_conc,
                                        self.ions['mw'].to_dict(),
                                        self.ions['valence'].to_dict())
        

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
        
        Cin_dict = self.ions.to_dict('index')
        u_Cin2 = {}
        val2 = {}
        MW2 = {}
        for c in Cin_dict.keys():
            u_Cin2[c] = Cin_dict[c]['units']
            val2[c] = Cin_dict[c]['valence']
            MW2[c] = Cin_dict[c]['mw']
        u_Cout2 = {}
        if units == None:
            u_Cout2 = u_Cin2
        else:
            for c in Cin_dict.keys():
                u_Cout2[c] = units
        
        u_Cin2 = {}
        for c in Cin_dict.keys():
            u_Cin2[c] = 'meq'
            
        temp_t = pd.Series(self.result.t)
        tmp_u = self.u_result[0, :, :]
        
        period_factor, u_in, u_out = conv_time('sec',period,'time', period)            
        idx = pd.Index(temp_t * period_factor, name = period)
        
        df_c = pd.DataFrame(tmp_u.T, index = idx, columns = self.names)
        
        tmp_df = df_c.copy(deep=True)

        C_out3, u_Cin2, u_Cout2 = conv_database(tmp_df, u_Cin2, u_Cout2, \
                                    conv_conc, MW2, val2)       
        
        col_name = C_out3.columns.tolist()
        
        new_col_name = []
        for n in col_name:
            new_col_name.append( n + ' ('+ u_Cout2[n] + '/L)')           

        C_out3.columns = new_col_name
        saved_name = output_file_name
        C_out3.to_excel(output_file_name, sheet_name = 'Cout', float_format='%.8f')        

        return saved_name
    
        
    def solve(self, t_eval=None):
        """ Returns (t, u)
        t = time values
        u = array of shape (phases, ions, time) 
        Liquid phase = 0 
        Resin phase = 1:
        """
        
       
        ### Alias some parameters for convenience ###

        Qm = np.float64(self.params['Qm']) # Resin Capacity by mass (meq/kg)
        RHOP = np.float64(self.params['RHOP']) # apparent resin density (kg/L)
        
        
        VR = np.float64(self.params['VR']) # resin phase volume (L)
        VL = np.float64(self.params['VL']) # liquid phase volume (L)
        
        nr = int(self.params['nr']) # number of radial collocation points
        Ds = np.float64(self.params['Ds']) # overall resin phase diffusion coefficient (cm**2/s)
        kL =  np.float64(self.params['kL'])   # overall film transfer coefficient (cm/s)
        rb = np.float64(self.params['rb']) # resin radius (cm)
        
        t_end = np.float64(self.params['t_end']) # experiment duration (seconds)
        
        C0 = self.ions['C0'].values
        
        
        ##########################
        ### DERIVED PARAMETERS ###
        ##########################
        

        
        # XXX: Might want some of these to be available after solving
        Q = RHOP * Qm # Resin capacity by volume 
        CT = C0.sum() # XXX: Need to modify for time varying feed concentrations
        
        NION = len(C0) 
        NEQ = (nr + 1) * (NION) # number of equations 
        
        # XXX: basically an enumeration. Probably should do this the Pythonic way
        LIQUID = 0
        RESIN = 1
        PRESAT = 0
        SURF = -1
        
        # divalent indexes
        valences = self.valences # XXX: Consider cleaning this up

        dv_idx = valences == 2
        ndv = len(valences[dv_idx])# number of divalent species
        
        # non chloride monovalent indexes
        mv_idx = valences == 1
        mv_idx[0] = False   # presaturant (reference ion) 
        
        # Equilibrium coefficients
        Kxc = self.ions['Kxc'].values  # separation factor vs presat, 
        # NOTE: the FIRST entry is the presaturant

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
        _, _, rootsr, Br, Wr = build_collocation(nr, 3)
        # Br is the radial, symmetric second derivative operator for a sphere
        # Wr is a vector of quadrature weights
        # rootsr are the coordinates of the collocation points

        # construct sparse block versions of operators
        Br_block = block_diag(*(NION*[Br])) 
        # that extra * in front unpacks the resulting list to block diag accepts the mats
##        self.Br_block = Br_block  # XXX: DEBUG
        Wr1_block = block_diag(*(NION*[Wr[:-1]]))
        
        
        ### Initialize concentrations ###
        u = np.zeros([(nr + 1), NION]) # equivalent concentration tensor
        u[LIQUID, :] = C0
        u[RESIN:, PRESAT] = Q  # resin initially loaded with PRESAT
        
        u0 = u.reshape(NEQ)  # initial concentration vector for solve_ivp call


        #########################
        ### Local Equilibrium ###
        #########################
        
        def calc_Ceq_dv(q_s, CT):
            # update Ceq with proper accounting for divalent ions
            Ceq = np.zeros(np.array([NION]))
            
            aa = (q_s[dv_idx] / Kxc[dv_idx]).sum() / q_s[PRESAT]**2
            bb = (q_s[mv_idx] / Kxc[mv_idx]).sum() / q_s[PRESAT] + 1
            cc = -CT
            Cc_eq = (-bb + np.sqrt(bb**2 - 4 * aa * cc))/(2*aa)
           
            Ceq[:] = q_s[:]/Kxc[:]*(Cc_eq/q_s[0])**valences[:]
                
            Ceq[PRESAT] = Cc_eq
            
            return Ceq
        
   
        def calc_Ceq_mv(q, CT):
            denoms = np.dot(alpha, q) # a vector of denominators
            Ceq = q * CT / denoms     # elementwise division 
            return Ceq
        
        if np.any(valences == 2):
            calc_Ceq = calc_Ceq_dv
            u[RESIN:, dv_idx] = Q/1e3 # so the quadratic formula doesn't explode
#           TODO: Test if alternate quadratic form can avoid this
            u[RESIN:, PRESAT] = Q - (ndv * Q / 1e3)  
            
        else:
            calc_Ceq = calc_Ceq_mv
            
        
        
        ###########################
        ### Solving Derivatives ###
        ###########################
        
        def diffun(T, u):
            """
            Calculate time derivatives at grid points
            """
            u = u.reshape([(nr + 1), NION])
        
            C = u[LIQUID, :]
            q = u[RESIN:, :]
            q_s = u[SURF, :]
            
            # update Ceq
            Ceq = calc_Ceq(q_s, CT)
            
            # Calculate flux terms
            J = - kL * (C - Ceq) # mass flux 
            a = 3 * VR / rb # surface area of a beads

            # Initialize arrays for derivatives
            du_dT = np.zeros([nr+1, NION])
            
            #### Liquid phase concentrations
            du_dT[LIQUID, :] = a * J / VL
            
            
            # diffusion in bead  
            q_swap = np.swapaxes(q, 0, 1).reshape(nr*NION)
            Br_q = Br_block.dot(q_swap)
            Br_q = Br_q.reshape([NION, nr])
            Br_q = np.swapaxes(Br_q, 0, 1)
            dq_dT =  Ds / rb**2 * Br_q   

            # intermediate term for dq_dT at bead surface        
            dq_dT_swap = np.swapaxes(dq_dT[:SURF, :],0,1).reshape(((nr-1)*NION))
            dq_dT_w = Wr1_block.dot(dq_dT_swap)
            dq_dT_w = dq_dT_w.reshape(NION)
            
            du_dT[RESIN:, :] = dq_dT
            du_dT[SURF, :] = (- J / rb - dq_dT_w) / Wr[-1]              

            du_dT = du_dT.reshape(NEQ) # solve_ivp requires this to be a vector
           
            return du_dT  
        
        
        #####################
        ### Solve Problem ###
        #####################
        
        T_final = t_end 
        Tvals = np.array([0, T_final]) # time points to solve
        
        T_eval = None
        if np.any(t_eval):  # specific times requested
            T_eval = t_eval 
        
        self.result = solve_ivp(diffun, Tvals, u0, method='BDF',
                                t_eval=T_eval, atol=1e-7, rtol=1e-5) 

        t = self.result.t 
        NT = len(t)
        u = self.result.y.reshape([(nr + 1), NION, NT])
        self.u_result = u
        
        ########################
        ### Check for Errors ###
        ########################
           
            
        if not np.allclose(u[1, :, :].sum(axis=0), Q, rtol=0.01):
            print('WARNING: Sum of resin concentrations is not Q!')
            
        return (t, u)
    

