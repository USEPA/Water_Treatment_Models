# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 13:56:50 2019

For more information on correlations for the film transfer coefficient see:
  Roberts, Paul V., Peter Cornel, and R. Scott Summers. 
  "External mass-transfer rate in fixed-bed adsorption." 
  Journal of Environmental Engineering 111, no. 6 (1985): 891-905.    
    
For more information on correlations for diffusivity in water see:    
  Li, Jianwei, and Peter W. Carr. 
  "Accuracy of empirical correlations for estimating diffusion coefficients in 
  aqueous organic mixtures." 
  Analytical chemistry 69, no. 13 (1997): 2530-2536.

@author: LHaupert
"""

import numpy as np


P_to_cP = 100 # viscosity unit conversion, poise to centipoise


# Dictionary of molecular volumes (mL/mol)
#chemspider.com
#ACD/PhysChem Suite
# Eventually we will want some kind of PFAS property database,
# but for now we'll use a dictionary.

MV_dict = {'PFBS':162.3,
           'PFOS':272.1,
           'PFBA':127.5,
           'PFHxS':217,
           'PFHxA':182.4,
           'PFHpA':210,
           'GenX':188.7,
           'PFOA':237.3,
           'PFNA':264.7,
           'PFDA':292.2} ## TODO: reference a file?


#  Water params
#  XXX:  Where did JB get these?

def viscosity(temp):
    ''' viscosity of water (g/cm-s) [AKA poise]
    temp in degC'''
#    return (1e-2)*(5.264789e-8*temp**4-1.207858e-5*temp**3+1.117294e-3*temp**2-\
#        5.661932e-2*temp+1.775586) 
    t = temp + 273.15
    return np.exp(-24.71 + (4209./t) + 0.04527 * t - (3.376e-5 * t**2))/100.
    

def density(temp):
    '''density of water gram/cm**3
    temp in degC'''  
    t = (temp + 273.15)/324.65
    return 0.98396*(-1.41768 + 8.97665*t - 12.2755 * t**2 + 7.45844 * t**3 - 1.73849 * t**4)
#    return 2.002888e-8*temp**3-6.334640e-6*temp**2+2.685668e-5*temp+1.000013 


def Hayduk_Laudie(MV, T):
    """
    Empirical correlation to estimate diffusion coefficient (cm2/s).
    T = temperature (C)
    MW = molecular weight of solvent (g/mol)
    MV = molar volume of solute (mL/mol)
    """
    
    mu = viscosity(T) * P_to_cP
    D_est = 13.26e-5 * (mu ** -1.14) * (MV ** -0.589)
    
    return D_est
   
    
def Wilke_Chang(MV, T,  MW=18.02, psi=2.6):
    """
    Empirical correlation to estimate diffusion coefficient (cm2/s).
    T = temperature (C)
    MW = molecular weight of solvent (g/mol)
    MV = molar volume of solute (mL/mol)
    psi = association constant (2.6 for water)
    """
    
    mu = viscosity(T) * P_to_cP
    kelvins = T + 273.15
    D_est = 7.4e-8 * np.sqrt(psi * MW) * kelvins * (MV ** -0.6) / mu
    
    return D_est
    

def simple_Gnielinski(rb, v, EBED, T, DL):
    """
    Returns estimate of film transfer coefficient (cm/s)
    rb = particle radius (cm)
    v = superficial flow velocity (cm/s)
    EBED = bed porosity
    T = temperature (deg. C)
    DL = contaminant diffusion coefficient in water (cm2/s)
    """
    
    rho_L = density(T)
    mu = viscosity(T)
    dP = 2* rb # bead diameter (cm)
    u = v / EBED  # interstitial linear flow velocity (cm/s)
    
    Re =  u * dP * (rho_L / mu)  # Reynolds number
    Sc = mu / rho_L / DL # Schmidt number
    Sh = (2 + 0.644 * Re**(1/2) * Sc**(1/3)) * (1 + 1.5 * (1- EBED)) # Sherwood number
    # print(Sh)
    
    # XXX: Make sure correlation is in valid range.
    if not (1 < Re < 100):
        print('WARNING: Reynolds number is outside the valid range for this correlation!')
    if not (0.6 < Sc < 1e4):
        print('WARNING: Schmidt number is outside the valid range for this correlation!')
    
    kL_est = Sh * DL / dP
    
    return kL_est


if __name__=='__main__':
    
    PFAS = 'PFOA'

    MV = MV_dict[PFAS] # molar volume (mL/mol)
    TEMPERATURE = 23 # deg C
    v = 0.272  # linear flow velocity (cm/s) [reported]
    EBED = 0.35 # bed porosity
    rb = 3.38e-2   # bead radius (cm)

    ### Estimate diffusion coefficient.
#    D_est = Wilke_Chang(MV, TEMPERATURE)
    D_est = Hayduk_Laudie(MV, TEMPERATURE)

    kL_est = simple_Gnielinski(rb, v, EBED, TEMPERATURE, D_est)

    print('Film Transfer Estimate: {:.2e} cm/s'.format(kL_est))


