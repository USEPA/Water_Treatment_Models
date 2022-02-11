import pandas as pd
import numpy as np

df = pd.read_csv('Chemicals.csv', header = 0)
#Load CSV File

def VapePressure(Chemical):
    
    ChemicalProperties = df[Chemical]
    #Column for Specific Chemical
    
    A = ChemicalProperties[0] #Get A Constant
    B = ChemicalProperties[1] #Get B Constant
    C = ChemicalProperties[2] #Get C Constant
    D = ChemicalProperties[3] #Get D Constant
    E = ChemicalProperties[4] #Get E Constant
    T = ChemicalProperties[5] #Get Temperature
    
    VP = np.exp(A + (B/T) + C*np.log(T) + D*T**E)
    #Formula from
    #Physical and Thermodynamic Properties of Pure Chemicals: Data Compilation
    
    return VP

print(VapePressure("Water"))