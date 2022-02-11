import pandas as pd
import numpy as np

df = pd.read_csv('Chemicals.csv', header = 0)
#Load CSV File

def LiquidDensity(Chemical):
   
    ChemicalProperties = df[Chemical]
    #Column for Specific Chemical
    
    A = ChemicalProperties[0] #Get A Constant
    B = ChemicalProperties[1] #Get B Constant
    C = ChemicalProperties[2] #Get C Constant
    D = ChemicalProperties[3] #Get D Constant
    T = ChemicalProperties[5] #Get Temperature
    
    if B*(1-(T/C)**D) != 0 and C !=0:
        #Prevent divide by 0
        
        LD = A/(B*(1-(T/C))**D)
        #Formula from 
        #Physical and Thermodynamic Properties of Pure Chemicals: Data Compilation
    
    return LD

print(LiquidDensity("Water"))