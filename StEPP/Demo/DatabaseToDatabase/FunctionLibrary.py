import pandas as pd
import numpy as np


def GetChemical(chemicals):
    
    #Read in Data
    ChemicalData = pd.read_csv('Chemicals.csv')
    
    
    #Define Lists
    index_list = []
    lst = []
    ChemicalList = ChemicalData['Chemical']
    
    
    #Convert ChemicalList from Object to List to make my life easier
    for i in ChemicalList:
        lst.append(i)
        
    
    #For element in chemicals, if element is in ChemicalList (lst), get index
    #of that element in ChemicalList
    for i in chemicals:
        for j in lst:
            if i==j:
                index_list.append(lst.index(i))
                
    
    #Get constants from database for chosen chemicals
    ChemicalVariableData = ChemicalData.loc[index_list, ['A','B','C','D','E']]
    
    
    #Create Dataframe of constants for later use
    ChemicalConstants = pd.DataFrame(ChemicalVariableData)
    
    
    #Return the value of those constants to put into formulas
    FormulaConstants = ChemicalVariableData.values
    
    
    return FormulaConstants


def liquiddensity(lst):
        
    for i in lst:
        calculated_ld = i/(i*(1-(t/i))**i)
        #Formula from 
        #Physical and Thermodynamic Properties of Pure Chemicals: Data Compilation
    
    return calculated_ld




def vaporpressure(lst):
    
    for i in lst:
        calculated_vp = np.exp(i + (i/t) + i*np.log(t) + i*t**i)
    #Formula from
    #Physical and Thermodynamic Properties of Pure Chemicals: Data Compilation
    
    return calculated_vp


def Results(chemicals, t=5):
    
    #Empty list to store calculations
    #List of Chemicals
    allproperties = []
    vaporpressures = []
    liquiddensities = []         
    chemicallist = chemicals        
    
    #Store variables to get calculations
    chemproperties = GetChemical(chemicals) 
    #Returns a list of chemical properties for each chemical ([[Chemial1_Properties],...,[ChemicalN_Properties]])
    
    #Add calculations to list
    for i in chemproperties:
        vaporpressures.append(vaporpressure(i))
        liquiddensities.append(liquiddensity(i))
    #For each chemical in list, add each vapor pressure to a list
    
    allproperties.append(vaporpressures)
    allproperties.append(liquiddensities)
    
    #Convert to datagrame
    properties = pd.DataFrame(data=allproperties, index=['Predicted Vapor Pressure', 'Predicted Liquid Density'], columns=chemicallist)
                                                                                 
    return properties


print(Results(['Water', 'O2','CO2']))







#print(ResultsDataframe(mychemicallist))

#properties = pd.DataFrame(data=GetChemical(['Water', 'O2']), index=['Water', 'O2'], columns=['A', 'B', 'C', 'D', 'E'])
#print(properties)

#water = GetChemical("Water")
#print(pd.DataFrame(data=liquiddensity(water), index=None, columns=None))

# def LiquidDiffusivity(Chemical):
    
#     mu = 0.00131 * 1000
#     #Viscosity of Water in g/cm/s to cp
    
#     Diffusivity = 13.26/((mu)**1.4*MolarVolume(Chemical)*100**2)
#     #Prediction of Diffusion Coefficients for Non-electrolytes in Dilute 
#     #Aqueous Solutions
    
#     return Diffusivity

