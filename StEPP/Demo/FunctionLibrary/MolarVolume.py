import pandas as pd
import numpy as np
import LiquidDensity

df = pd.read_csv('Chemicals.csv', header = 0)
#Load CSV File

print(LiquidDensity.LiquidDensity("Water"))

def MolarVolume(Chemical):
    
    ChemicalProperties = df[Chemical]
    #Column for Specific Chemical
    
#    ChemicalMass = ChemicalProperties[0]
#    LiqDen = LD(Chemical)
    
#    Volume = ChemicalMass/LiquidDensity

#    return Volume

print(MolarVolume("Water"))
