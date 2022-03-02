import pandas as pd
import numpy as np

t=5

def GetChemical(chemical):
    
    #Read Data
    ChemicalData = pd.read_csv('Chemicals.csv')
    
    #Get Constants
    #Get Index of Row we want
    lst = []
    ChemicalObject = ChemicalData['Chemical']
    ChemicalList = ChemicalObject.values
    
    for i in ChemicalList:
        lst.append(i)
    
    indexer = lst.index(chemical)
    
    ChemicalVariableData = ChemicalData.loc[indexer, ['A','B','C','D','E']]
    
    #Create Dataframe
    #Put constants in dataframe_constants
    
    ChemicalConstants = pd.DataFrame(ChemicalVariableData)
    
    #Assign data for function
    FormulaConstants = ChemicalVariableData.values
    
    return FormulaConstants



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