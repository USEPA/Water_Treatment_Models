import pandas as pd
import numpy as np

def LiquidDensity(lst):
        
        for i in lst:
            LD = i/(i*(1-(t/i))**i)
        #Formula from 
        #Physical and Thermodynamic Properties of Pure Chemicals: Data Compilation
    
    return LD




def VapePressure(lst):
    
    for i in lst:
        VP = np.exp(i + (i/t) + i*np.log(t) + i*t**i)
    #Formula from
    #Physical and Thermodynamic Properties of Pure Chemicals: Data Compilation
    
    return VP



# def MolarVolume(Chemical):
    
#     ChemicalProperties = df[Chemical]
#     #Column for Specific Chemical
    
#     ChemicalMass = ChemicalProperties[0]
#     LD = LiquidDensity.LiquidDensity(Chemical)
    
#     Volume = ChemicalMass/LiquidDensity
#     #Formula from 
#     #Physical and Thermodynamic Properties of Pure Chemicals: Data Compilation

#     return Volume

# print(MolarVolume("Water"))



# df = pd.read_csv('Chemicals.csv', header = 0)
# #Load CSV File

# def LiquidDiffusivity(Chemical):
    
#     mu = 0.00131 * 1000
#     #Viscosity of Water in g/cm/s to cp
    
#     Diffusivity = 13.26/((mu)**1.4*MolarVolume(Chemical)*100**2)
#     #Prediction of Diffusion Coefficients for Non-electrolytes in Dilute 
#     #Aqueous Solutions
    
#     return Diffusivity