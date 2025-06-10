'''
PAC PSDM Model

This is a CFPSDM model under development for modeling powdered activated carbon.

This code is still being tested.

Authors:
Jonathan Burkhardt
Levi Haupert
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root

import sys
import os
cwd = os.getcwd()
split_cwd = cwd.split("\\")
sys_path = "/".join(split_cwd[:-1]+["PSDM/psdm"])
sys.path.append(sys_path)

from PSDM_functions import calc_solver_matrix, viscosity, density

try:
    import mkl
    mkl.set_num_threads(1)
except:
    pass


#### unit conversion helpers
lpg = 3.785411784 # liter per gallon conversion
cm_per_ft = 2.54 * 12
min_per_day = 24 * 60
gm_per_lb = 453.59237


### time conversions might need to end up in seconds
time_convert = {'min': 60, 'mins': 60, 
                'hr': 60*60, 'hrs': 60*60, 
                'day': min_per_day*60, 'days': min_per_day*60,
                's': 1, 'sec': 1, 'secs':1
                } ## converts to seconds

length_convert = {'cm': 1, 'm': 100, 'mm': 0.1,
                  'in': 2.54, 'ft': cm_per_ft,
                  }  ## converts to cm

volume_convert = {'ml': 1, 'cm3': 1,
                  'l': 1e3, 
                  'gal': lpg * 1e3,
                  'ft3': length_convert['ft']**3,
                  'm3': length_convert['m']**3}

flow_convert = {'ml/min': 1/time_convert['min'], 'ml/s': 1,
                'gpm': lpg * 1e3/time_convert['min'], 'gal/min': lpg * 1e3/time_convert['min'],
                'lpm': 1e3/time_convert['min'], 'l/min': 1e3/time_convert['min'],
                'mgd': 1e6 * lpg * 1e3 / time_convert['day'],
                'ft3/s': length_convert['ft']**3/time_convert['s'], 
                'cfs': length_convert['ft']**3/time_convert['s'],
                'm3/s': length_convert['m']**3/time_convert['s']
                } ## convert to ml/min

conc_convert = {'ug/L': 1, 'ug': 1,
                'ng/L': 1e-3, 'ng': 1e-3,
                'mg/L': 1e3, 'mg': 1e3
                } ## convert to ug/L

# conc_convert

class PAC_CFPSDM():
    def __init__(self, contactor_df, pac_df, compounds_df, help_print=False, **kw):

        ## force correct unit types, primarily used from R Shiny app
        for idx in contactor_df.index:
            if idx in ['format']:
                contactor_df.loc[idx, 'value'] = str(contactor_df.loc[idx, 'value'])
            else:
                ## pretty much everything should tolerate being a float
                contactor_df.loc[idx, 'value'] = float(contactor_df.loc[idx, 'value'])#.astype('float64')

        for idx in pac_df.index:
            ## everything should be a float
            pac_df.loc[idx, 'value'] = float(pac_df.loc[idx, 'value'])

        for idx in compounds_df.index:
            if idx != 'C0_units':
                ## C0_units should remain a string
                compounds_df.loc[idx] = compounds_df.loc[idx].astype('float64')
        ### end unit correction for R
        
        

        self.help_print = help_print
        self.errors = 0 ## count of errors

        self.contactor_index = [i.lower() for i in contactor_df.index]
        contactor_df.index = self.contactor_index   ## resets the index column to all lower case
        self.pac_index = [i.lower() for i in pac_df.index]
        pac_df.index = self.pac_index               ## resets the index column to all lower case

        ## initiate collocation
        self.nc = kw.get('nr', 5)  #set number of radial points
        self.mc = kw.get('nz', 8) #set number of axial points, or 12, not needed for PAC, but needed for calc_solver_matrix
        self.nz = self.mc * 1
        solver_data = calc_solver_matrix(self.nc, self.mc, 1)
        self.wr = solver_data['wr']
        self.az = solver_data['az']
        self.br = solver_data['br']
        if self.mc != solver_data['mc']:
            ''' corrects for OCFE change to nz'''
            self.mc = solver_data['mc']
        self.nd = self.nc - 1

        #set up temperature dependant values
        if 'temperature' in self.contactor_index:
            self.temp = contactor_df.loc['temperature', 'value'] 
        else:
            self.temp = 20 ### assumes 20 degC if no temperature is provided
        self.vw = viscosity(self.temp)      ### gm/cm-s
        self.dw = density(self.temp)        ### gm/cm**3
                
        self.time_type = kw.get('time_type', 'min')
        self.t_mult = time_convert[self.time_type]

        
        ### PAC information #######################################################################################
        if 'density' in pac_df.index:
            ## expected gm/mL
            self.density = pac_df.loc['density','value'] ### may need to add unit conversion
        else:
            self.density = np.nan
            self.errors += 1
            print('Please provide "density" in pac_df')

        if 'porosity' in pac_df.index:
            ## unitless
            self.porosity = pac_df.loc['porosity', 'value']
        else:
            self.porosity = np.nan
            self.errors += 1
            print('Please provide "porosity" in pac_df')

        if 'radius' in pac_df.index:
            ## converted to cm
            self.pac_radius = pac_df.loc['radius', 'value'] * length_convert[pac_df.loc['radius','units'].lower()]
        else:
            self.pac_radius = np.nan
            self.errors += 1
            print('Please provide "radius" in pac_df')
        

        ## column sizing and characteristics ############################################################################
        if 'length/diameter' in self.contactor_index:
            ## converted to cm
            self.length = contactor_df.loc['length/diameter', 'value'] * length_convert[contactor_df.loc['length/diameter', 'units'].lower()]
        else:
            print('Please provide "length/diameter" in contactor_df')
            self.errors += 1

        self.format = 'square' ### assumes square by default
        if 'format' in self.contactor_index:
            self.format = contactor_df.loc['format', 'value']

        self.width = np.nan  ## may not be needed for circular, default nan
        if self.format == 'square':
            self.width = self.length
        elif self.format == 'rectangular':
            if 'aspect ratio' in self.contactor_index:
                self.aspect_ratio = contactor_df.loc['aspect ratio', 'value'] ## L/W
                self.width = self.length / self.aspect_ratio
            else:
                self.aspect_ratio = 1
                self.width = self.length / self.aspect_ratio
                print('Warning: No "aspect ratio" was provided, assuming 1')

        if 'height' in self.contactor_index:
            self.height = contactor_df.loc['height', 'value'] * length_convert[contactor_df.loc['height', 'units'].lower()]
        else:
            print('Please provide "height" in contactor_df')
            self.height = np.nan
            self.errors += 1

        self.volume = np.nan
        if 'volume' in self.contactor_index:
            ## check if other values were already provided
            
            if self.format == 'circular':
                test = np.isnan(self.length * self.height)
            else:
                test = np.isnan(self.length * self.width * self.height)

            if test:
                ## Only worries about volume if other values aren't provided, otherwise lets it get calculated.
                self.volume = contactor_df.loc['volume', 'value'] * volume_convert[contactor_df.loc['volume', 'units'].lower()]

                ## use a specified volume to calculate other terms
                if np.isnan(self.height):
                    if ~np.isnan(self.length * self.width):
                        ## should capture square or rectangular
                        self.errors -= 1
                        self.height = self.volume / (self.length * self.width)
                    elif self.format == 'circular' and ~np.isnan(self.length):
                        self.errors -= 1
                        self.height = self.volume / (np.pi/4 * (self.length ** 2))
                if np.isnan(self.length) and ~np.isnan(self.height):
                    if self.format == 'circular':
                        self.errors -= 1
                        self.length = np.sqrt(self.volume * 4 / (np.pi * self.height))
                    else:
                        self.format == 'square' ## forces square
                        self.errors -= 1
                        self.length = np.sqrt(self.volume/self.height)
                        self.width = self.length * 1

        self.flow = np.nan
        if 'flow' in self.contactor_index:
            self.flow = contactor_df.loc['flow', 'value'] * flow_convert[contactor_df.loc['flow', 'units'].lower()]
        else:
            self.error += 1
            print('Please provide "flow" in contactor_df')

        if np.isnan(self.volume):
            if self.format == 'circular':
                self.volume = self.height * np.pi/4 * (self.length**2)
            else:
                ## length and width should be already provided
                self.volume = self.length * self.width * self.height

        if 'hrt' in self.contactor_index:
            self.hrt = contactor_df.loc['hrt', 'value'] * time_convert[contactor_df.loc['hrt', 'units'].lower()]
        else:
            self.errors += 1
            self.hrt = np.nan
            print('Please provide "HRT" for hydraulic residence time in contactor_df')

        if 'crt' in self.contactor_index:
            self.crt = contactor_df.loc['crt', 'value'] * time_convert[contactor_df.loc['crt', 'units'].lower()]
        else:
            self.errors += 1
            self.crt = np.nan
            print('Please provide "CRT" for carbon residence time in contactor_df')

      
        self.duration = kw.get('duration', np.max([self.hrt, self.crt])) ## set duration to the max of HRT and CRT unless provided, ## units in seconds?
        self.time_incr = time_convert['sec']

        if 'pac dosage' in self.contactor_index:
            self.dosage = contactor_df.loc['pac dosage', 'value']  ## need to convert? or assume mg/L always
        else:
            self.errors += 1
            self.dosage = np.nan
            print('Please provide "PAC Dosage" in contactor_df')

        self.TIS = np.nan
        if 'tanks in series' in self.contactor_index:
            if contactor_df.loc['tanks in series', 'value']:
                self.TIS = contactor_df.loc['tanks in series', 'units']  ## should be a number in this case


        self.compounds_df = compounds_df.copy()
        self.convert_array = np.array([conc_convert[i] for i in self.compounds_df.loc['C0_units']])
        
        ### check on correlations for kf, Ds, Dp   ### TODO!!

        self._update_values()  ### using a function, so I can set up an easier iterator on self.dosage, like the GUI has.
        

        self.ncomp = len(self.compounds_df.columns)
        self.effluent_locator = [(self.nc+1)*(i+1)-1 for i in range(self.ncomp)]
        ### address variable influent?? Maybe next step.


        ### End matter ############################################################################
        
        if self.help_print:
            print(f'PAC density: {self.density:.3f} g/mL\tPorosity: {self.porosity:.3f}\tRadius: {self.pac_radius:.3e} cm')
            print(f'Temp: {self.temp} degC - Viscosity: {self.vw:.3e} g/cm-s \tDensity: {self.dw:.3f} gm/cm**3')
            print(f'Time multiplier: {self.t_mult}')
            print(f'{self.format} - Length {self.length:.2f} cm, Width {self.width:.2f} cm -- Height: {self.height:.2f} cm')
            print(f'Volume: {self.volume:,.2f} mL')
            print(f'Flowrate: {self.flow:,.2f} mL/min')
            print(f'Hydraulic Residence Time (HRT): {self.hrt:.2f} min')
            print(f'Carbon Residence Time (HRT): {self.crt:.2f} min')
            print(f'PAC Dosage: {self.dosage:.2f} mg/L')

            # print(self.compounds_df)

        if self.errors > 0:
            print(f'{self.errors} error(s) was/were found in the input, please review')

    def _update_values(self):
        self.epsilon = 1 - self.dosage/(self.dw * 1e6)

        self.compounds_df.loc['qe'] = self.compounds_df.loc['K'] * (self.compounds_df.loc['C0'] * self.convert_array)**self.compounds_df.loc['1/n'] ## should be in ug/g

        self.compounds_df.loc['Dgs'] = self.density * 1e3 * self.compounds_df.loc['qe'] * (1 - self.epsilon) / (self.epsilon * (self.compounds_df.loc['C0'] * self.convert_array))
        self.compounds_df.loc['Dgp'] = self.porosity * (1 - self.epsilon) / (self.epsilon)
        self.compounds_df.loc['Dg'] = self.compounds_df.loc['Dgs'] + self.compounds_df.loc['Dgp']

        self.compounds_df.loc['St'] = self.compounds_df.loc['kf'] * self.hrt * self.porosity * (1 - self.epsilon) / (self.epsilon * self.pac_radius) #  
        self.compounds_df.loc['Bi'] = self.pac_radius * self.compounds_df.loc['kf'] * (1 - self.epsilon) / (self.epsilon * (self.compounds_df.loc['Ds'] * self.compounds_df.loc['Dg'] + self.compounds_df.loc['Dp'] * self.compounds_df.loc['Dgp']))

        self.compounds_df.loc['Z'] = self.compounds_df.loc['Dgp'] * (self.compounds_df.loc['Dp'] - self.compounds_df.loc['Ds']) * self.hrt / (self.pac_radius * self.pac_radius * self.compounds_df.loc['Dg']) # * 60
        self.compounds_df.loc['X'] = self.compounds_df.loc['Ds'] * self.hrt / (self.pac_radius * self.pac_radius)  # * 60 ## unitless

        self.compounds_df.loc['Ye'] = self.compounds_df.loc['qe'] + self.porosity / (self.density * 1e3) * self.compounds_df.loc['C0'] * self.convert_array  ## ug/g
        
        if self.help_print:
            print(self.compounds_df)
            print('Epsilon', self.epsilon)
    
    def get_analytical_sol(self, t):
        ## Film Transfer only assumption
        a = 1 / (self.epsilon * self.hrt) + (3 * self.compounds_df.loc['kf'].values.astype('float64') * (1 - self.epsilon))/(self.pac_radius * self.epsilon)

        b = self.compounds_df.loc['C0'].values.astype('float64') / (self.epsilon * self.hrt) 

        c = a * self.compounds_df.loc['C0'].values.astype('float64') - b
        d = np.exp(np.outer(a, t))

        ceff = (np.divide(c, d.T) + b) / a

        return ceff
    
    def get_analytical_sol2(self):
        def find_qn(n, A):
            """ WARNING: only works for positive values of A """
            roots = []
            guess = np.pi - 0.000001 

            def f(x):
                return np.tan(x) - (3*x / (3+A * x**2))

            for N in range(n):
                if roots:
                    guess += np.pi  # move to next interval
                roots.append(brentq(f, guess, guess + np.pi/2)) # XXX: Will fail if A is negative
            return roots
        
        def find_pn(n, A, B):
            """ WARNING: only works for positive values of A """
            roots = []
            guess = np.pi / 2 + 0.000001 

            def f(x):
                numer = 3 - (1-A)/A * x**2 / B
                denom = 3 + (1-A)/A * (B-1)/B * x**2
                return numer/denom - np.tan(x)/x 

            for N in range(n):
                
                if not roots:  # XXX: Try some hackish things to see if the first root is split.
                
                    try:
                        roots.append(brentq(f, guess, (guess + np.pi*0.9999)/2))
                    except:
                        # print('Root failure 1')
                        pass
                    
                    try:
                        roots.append(brentq(f, (guess + np.pi*0.9999)/2, guess + np.pi*0.9999))     
                    except:
                        # print('Root failure 2')
                        pass
                
                    
                if roots:
                    guess += np.pi  # move to next interval
                roots.append(brentq(f, guess, guess + np.pi*0.9999)) 
                
           
            return roots

        def finiteSphere(t, n, D, ra, A):
            """Return a function for finding M$_T$ / M$_eq$ for a sphere immersed in a
            well stirred bath of finite size.
            See Crank eq 6.30
            """

            q = find_qn(n, A)
            q = np.array(q)
            # Making a 2D array here should be faster than Sympy expamsion 
            q = q[None, :]
            t = t[:, None]
        
            B = 6 * A * (A + 1)
            num = B * np.e ** (-D * t * q**2 / ra**2)
            denom = 9 + 9*A + q**2 * A**2
            
            sum_term = (num/denom).sum(axis=1)
            
            f = 1 - sum_term
            
            return f
        
        def helfferich_film(t, rb, kL, Vw, Vx, kAxw):
            """
            Return fractional progress to equilibrium over time
            A is contaminant
            B is presaturant
            [I'm sorry, this seems backwards notationally to what is typically used]
            """

            #    kAxw = qAf / cAf    # Distribution coefficient
            #    kAxw = KAB * Q / cBi  # Trace contaminant approximation


            LAM = kAxw * Vx / Vw
        
            # prefactor
            ex_fact = -3 * kL / rb / kAxw   # eq 30, 31, except for t

            f = 1 - np.exp((LAM+1)*(ex_fact * t)) # eq 67
        
            return f
        
        def finiteHSDM(t, kf, Ds, rb, Vw, Vx, Kd, nterms=5):
            """Return a function for finding M$_T$ / M$_eq$ for a sphere immersed in a
            well stirred bath of finite size.
            See Perry's Handbook chapter 16
            """
            n = nterms
            Bi = kf * rb / Ds / Kd
            LAM = 1 / (1 + (Vw / Vx / Kd)) 
            
            q = find_pn(n, LAM, Bi)
            q = np.array(q)
            # Making a 2D array here should be faster than Sympy expansion 
            q = q[None, :]
            t = t[:, None]
        
            num = 6 * np.e ** (- Ds * t * q**2 / rb**2)
            denom = 9 * LAM / (1-LAM) + (1-LAM)*q**2 -(5*LAM + 1) * q**2 / Bi + (1-LAM) * q**4 / Bi**2
            
            sum_term = (num/denom).sum(axis=1)
            
            f = 1 - sum_term
            
            return f

        self.duration_old = self.duration * 1
        self.duration = 3000 * 60 ## make sure it has reached equilibrium

        results = self.run_PAC_PSDM()

        # Otherwise we are going to have to trust the endpoint of the numerical sim
        times = self.y.t

        Ye = 3 * self.wr.dot(self.y.y.reshape(self.nc+1, self.ncomp, len(times))[:-1, :, -1])
        Ce = self.y.y.reshape(self.nc+1, self.ncomp, len(times))[-1, :, -1]
        # qe_s = qe(Ce, K, xn)

        Kd = Ye / Ce * self.density
        # print('\nYe', Ye, '\nCe', Ce,'\nKd', Kd,)# '\nqe', qe_s)
        # print(self.compounds_df.loc['Ye'].values)


        ## Crank
        Deff = self.compounds_df.loc['Ds'].values + self.porosity * self.compounds_df.loc['Dp'].values / (self.density * 1e3 * Kd)
        
        Bi = self.pac_radius * self.compounds_df.loc['kf'].values / Deff / Kd

        A = 1 / (self.dosage / (self.density * 1e6)) / Kd
        depletion = (1 / (1 + A)).reshape((self.ncomp, 1))

        f = np.array([finiteSphere(times, 500, Deff[i], self.pac_radius, A[i]) for i in range(self.ncomp)])
        F = 1 - depletion * f 

        ## Helfferich
        VL = 1 ## L
        Vp = self.dosage / (self.density * 1e6)
        f = np.array([helfferich_film(times, self.pac_radius, self.compounds_df.loc['kf'].values[i], VL, Vp, Kd[i]) for i in range(self.ncomp)])
        F_h = 1 - f * depletion

        ### HSDM
        hsdm_error = False
        Deff = self.compounds_df.loc['Ds'].values
        f = []
        for i in range(self.ncomp): 
            try:
                f_part = list(finiteHSDM(times, self.compounds_df.loc['kf'].values[i], Deff[i], self.pac_radius, VL, Vp, Kd[i], nterms=10))
                f.append(f_part)
            except:
                hsdm_error = True
                
        if not hsdm_error:
            f = np.array(f) 
            F_hsdm = 1 - f * depletion 
        else:
            F_hsdm = np.zeros((self.ncomp, len(times))) * np.nan

        ## reset values
        self.duration = self.duration_old * 1

        return times, F, F_h, F_hsdm



    def run_PAC_PSDM(self):

        self.y = np.zeros((self.nc+1, self.ncomp))
        self.yprime = np.zeros((self.nc + 1, self.ncomp))
        self.y[-1] = 1 ### set liquid balance

        self.cPore = np.zeros((self.nc, self.ncomp))

        ## create some convenience arrays
        xn_array = self.compounds_df.loc['1/n'].values.astype('float64')
        ye_array = self.compounds_df.loc['Ye'].values.astype('float64')
        x_array = self.compounds_df.loc['X'].values.astype('float64')
        z_array = self.compounds_df.loc['Z'].values.astype('float64')
        bi_array = self.compounds_df.loc['Bi'].values.astype('float64')
        cin_array = self.compounds_df.loc['C0'].values.astype('float64') * self.convert_array ## convert to ug/L
        st_array = self.compounds_df.loc['St'].values.astype('float64')

        helper_array = np.ones((self.nc, self.ncomp))

        RHO = self.density * 1000 
        Kf = self.compounds_df.loc['K'].values.astype('float64')
        n = 1 / self.compounds_df.loc['1/n'].values.astype('float64')
        kL = self.compounds_df.loc['kf'].values.astype('float64')
        E = self.epsilon
        Ep = self.porosity
        Dp = self.compounds_df.loc['Dp'].values.astype('float64')
        Ds = self.compounds_df.loc['Ds'].values.astype('float64')
        B = self.br
        ra = self.pac_radius
        W = self.wr
        a_s = 3 / ra

        self.flag = 'PSDM'
        def diffun(t, y0):
            # print(f'{t/self.duration:.3f} %')
            y = y0.copy().reshape((self.nc + 1, self.ncomp))
            yprime = np.zeros((self.nc + 1, self.ncomp))            


            C = y[-1]
            Y = y[:-1]

            # calc local equilibrum pore conc given for given YY
            # NOTE: assume Y / RHO ~= q
            q = Y / RHO / Kf
            q[q < 0.0] = 0.0 # avoid taking a negative number to the n power
            Cp = q ** n 
            
            # SURFACE FLUX TERM
            J = - kL * (C - Cp[-1])    
            
            # Liquid phase
            dC_dt = (1 - E) / E * J * a_s 
        
            # diffusion in bead
            dY_dt = (Ep * (Dp - Ds) * B.dot(Cp) + Ds * B.dot(Y)) / ra**2
            
            # Intermediate term for boundary at bead surface
            dY_dt_w = W[:-1].dot(dY_dt[:-1])
            
            # Boundary condition at bead surface
            dY_dt[-1] = -(J/ra + dY_dt_w) / W[-1]    
        
            # output derivatives
            # du_dt = np.zeros(y.shape)
            # du_dt[-1] = dC_dt
            # du_dt[:-1] = dY_dt

            du_dt = np.append(dY_dt, dC_dt.reshape((1, self.ncomp)), axis=0)

            yprime = du_dt

            return yprime.flatten()


        self.y = solve_ivp(diffun, 
                      (0, self.duration),
                      self.y.flatten(),
                      solver='BDF',
                      )
        
        times = self.y.t / time_convert['min'] ### converts times to minutes

        effluent = np.multiply(self.compounds_df.loc['C0'].values.astype('float64'), self.y.y.reshape(self.nc+1, self.ncomp, len(times))[-1].T)
        out_df = pd.DataFrame(effluent, columns=self.compounds_df.columns, index=times)

        return out_df
    
    
    def run_multi_dosage(self, dosages):

        ## assumes dosages is an iterable
        data_dict = {i: [] for i in dosages}

        ## create stored values
        self.orig_dosage = self.dosage * 1
        self.orig_duration = self.duration * 1

        self.duration = 1000 * 60 ## duration in seconds

        for dose in dosages:
            self.dosage = dose * 1
            self._update_values()

            data_dict[dose] = self.run_PAC_PSDM()


        ### reset values
        self.dosage = self.orig_dosage * 1 ## reset to original information
        self.duration = self.orig_duration * 1

        return data_dict
        # return pd.DataFrame(data_dict) # Test -CDS
    
    def multi_dosage_analyzer(self, data_dict, target_HRT):
        '''
        returns dictionary of dataframes with results
        '''
        ## data_dict expected to have keys of dosage, and data for concentration lost by HRT
        ### Do we just have this analyzer run the function by itself? 

        if type(target_HRT) == int or type(target_HRT) == float:
            target_HRT_array = np.array([target_HRT])
        else:
            target_HRT_array = np.array(target_HRT) ## makes this an array, assumes array already or list
        
        out_dict = {self.compounds_df.columns[i]: pd.DataFrame(index=data_dict.keys(), columns=target_HRT_array) for i in range(self.ncomp)}

        for key in data_dict.keys():
            y_df = data_dict[key]
            
            conv_t = y_df.index ## in minutes by default
            
            for i in range(self.ncomp):
                compound = y_df.columns[i]
                ## calculate the concentration at a given HRT for each compound
                concs_calced = np.interp(target_HRT_array, conv_t, y_df[compound]) * self.convert_array[i]

                out_dict[compound].loc[key] = concs_calced * 1

        return out_dict

    def _run_multi_doseR(self, dosages, target_HRT):
        if type(target_HRT) == int or type(target_HRT) == float:
            target_HRT_array = np.array([target_HRT])
        else:
            target_HRT_array = np.array(target_HRT) ## makes this an array, assumes array already or list

        ## assumes dosages is an iterable
        data_dict = {float(i): [] for i in dosages}

        ## create stored values
        self.orig_dosage = self.dosage * 1
        self.orig_duration = self.duration * 1

        self.duration = 1000 * 60 ## duration in seconds

        for dose in dosages:
            self.dosage = dose * 1
            self._update_values()

            data_dict[dose] = self.run_PAC_PSDM()

        out_dict = {self.compounds_df.columns[i]: pd.DataFrame(index=data_dict.keys(), columns=target_HRT_array) for i in range(self.ncomp)}

        for key in data_dict.keys():
            y_df = data_dict[key]
            
            conv_t = y_df.index ## in minutes by default
            
            for i in range(self.ncomp):
                compound = y_df.columns[i]
                ## calculate the concentration at a given HRT for each compound
                concs_calced = np.interp(target_HRT_array, conv_t, y_df[compound]) * self.convert_array[i]

                out_dict[compound].loc[key] = concs_calced * 1

        
        # out_df = pd.DataFrame(index=pd.MultiIndex.from_tuples((i, j) for i in y_df.columns for j in dosages), columns=target_HRT_array)
        # for i in y_df.columns:
        #     for j in dosages:
        #         for k in target_HRT_array:
        #             out_df.loc[(i, j), k] = out_dict[i].loc[j, k] * 1
       
        ### reset values
        self.dosage = self.orig_dosage * 1 ## reset to original information
        self.duration = self.orig_duration * 1

        return out_dict


    def HRT_calculator_for_dosage(self, target_conc: float, target_units, dosage_trials=np.arange(5, 151, 20), influent_c0=np.arange(5, 26, 5), conc_units='same') -> dict:
        '''
        returns dictionary of dataframes with results of HRT for a given dosage
        '''
        ## save original data
        original_data = self.compounds_df.copy()
        original_duration = self.duration * 1
        original_conv_array = self.convert_array * 1

        self.duration = 3000 * 60

        target_conc_update = target_conc * conc_convert[target_units.lower()]
        
        if conc_units.lower() == 'same':
            ## sets conc_units to be the same as target_units
            conc_multiplier = np.ones(self.ncomp) * conc_convert[target_units.lower()]
            conc_units_update = [target_units.lower()] * self.ncomp ## make it the same size as the compounds_df shape
            self.convert_array = np.array([conc_convert[i] for i in conc_units_update])
        else:
            ### assumes units
            # 'ng', 'ug', 'mg'
            conc_multiplier = np.ones(self.ncomp) * conc_convert[conc_units.lower()]
            conc_units_update = [conc_units.lower()] * self.ncomp ## make it the same size as the compounds_df shape
            self.convert_array = np.array([conc_convert[i] for i in conc_units_update])


        out_dict = {self.compounds_df.columns[i]: pd.DataFrame(index=dosage_trials, columns=influent_c0) for i in range(self.ncomp)}

        self.compounds_df.loc['C0_units'] = conc_units_update

        for c0 in influent_c0: ### needs this to be an iterable
            for dosage in dosage_trials:

                ## change the inputs
                self.compounds_df.loc['C0'] = c0 * 1.
                self.dosage = dosage
                self._update_values()

                
                data = self.run_PAC_PSDM()
   
                for i in range(self.ncomp):
                    comp = data.columns[i]
                    sub_data = data[comp].values * self.convert_array[i]
                    sub_data_filtered = sub_data[np.where(sub_data <= target_conc_update)]

                    if len(sub_data_filtered) > 0:
                        sol = np.interp(target_conc_update, sub_data.flatten()[::-1], data.index[::-1])

                        out_dict[comp].loc[dosage, c0] = sol * 1.
                
        # plt.plot([0, self.duration/60], np.ones(2) * target_conc_update[0], lw=4, ls='--')

        ## reset compounds_df
        self.compounds_df = original_data.copy()
        self.duration = original_duration * 1
        self.convert_array = original_conv_array * 1

        return out_dict
        

    def _R_HRT_calculator_for_dosage(self, target_conc: float, target_units='ng', dosage_trials=np.arange(5, 150 + 1, 20), conc_units='same') -> pd.DataFrame:
        '''
        returns dataframe with results of HRT for a given dosage
        '''
        ## save original data
        original_data = self.compounds_df.copy()
        original_duration = self.duration * 1
        original_conv_array = self.convert_array * 1

        self.duration = 3000 * 60

        target_conc_update = target_conc * conc_convert[target_units.lower()]
        # print(target_conc_update)
        
        if conc_units.lower() == 'same':
            ## sets conc_units to be the same as target_units
            pass
            # conc_multiplier = np.ones(self.ncomp) * conc_convert[target_units.lower()]
            conc_units_update = [target_units.lower()] * self.ncomp ## make it the same size as the compounds_df shape
            # self.convert_array = np.array([conc_convert[i] for i in conc_units_update])
        else:
            ### assumes units
            # 'ng', 'ug', 'mg'
            conc_multiplier = np.ones(self.ncomp) * conc_convert[conc_units.lower()]
            conc_units_update = [conc_units.lower()] * self.ncomp ## make it the same size as the compounds_df shape
            self.convert_array = np.array([conc_convert[i] for i in conc_units_update])

         ### adjust target conc into appropriate unit

        ### calculates 
        out_df = pd.DataFrame(columns=self.compounds_df.columns, index=dosage_trials)

        self.compounds_df.loc['C0_units'] = conc_units_update

        # for c0 in influent_c0: ### needs this to be an iterable
        for dosage in dosage_trials:

            ## change the inputs
            # self.compounds_df.loc['C0'] = c0 * 1.
            self.dosage = dosage
            self._update_values()

            data = self.run_PAC_PSDM()

            for i in range(self.ncomp):
                comp = data.columns[i]
                sub_data = data[comp].values * self.convert_array[i]
                sub_data_filtered = sub_data[np.where(sub_data <= target_conc_update)]
                
                # plt.plot(data.index, data[comp].values * conc_multiplier[i], label=f'{c0} - {dosage}')

                if len(sub_data_filtered) > 0:
                    sol = np.interp(target_conc_update, sub_data.flatten()[::-1], data.index[::-1])

                    out_df.loc[dosage, comp] = sol * 1.
                

        ## reset compounds_df
        self.compounds_df = original_data.copy()
        self.duration = original_duration * 1
        self.convert_array = original_conv_array * 1

        return out_df
    
        
def R_run_PAC(contactor_df, pac_df, compounds_df, nr):
    ### move initial column to index for python version
    contactor_df.index = contactor_df['name']
    contactor_df = contactor_df[['value', 'units']]
    
    pac_df.index = pac_df['name']
    pac_df = pac_df[['value', 'units']]

    compounds_df.index = compounds_df[compounds_df.columns[0]]
    compounds_df = compounds_df[compounds_df.columns[1:]]

    for idx in compounds_df.index:
        if idx == 'C0_units':
            pass
            # print(compounds_df.loc[idx])
            # print('yup')
    # print(compounds_df)
   
    pac_mod = PAC_CFPSDM(contactor_df, pac_df, compounds_df, nr=nr)
    # print(pac_mod.run_PAC_PSDM())

    return pac_mod.run_PAC_PSDM() * pac_mod.convert_array


