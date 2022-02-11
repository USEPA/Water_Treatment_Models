import pandas as pd
import numpy as np
import MolarVolume

df = pd.read_csv('Chemicals.csv', header = 0)
#Load CSV File

def LiquidDiffusivity(Chemical):
    
    mu = 0.00131 * 1000
    #Viscosity of Water in g/cm/s to cp
    
    Diffusivity = 13.26/((mu)**1.4*MolarVolume(Chemical)*100**2)
    #Prediction of Diffusion Coefficients for Non-electrolytes in Dilute 
    #Aqueous Solutions
    
    return Diffusivity


