import pandas as pd
import numpy as np


#-----------------------------------------------------------------------------#
#This functions purpose is to obtain the chemical properties for a number of
#n chemicals. These properties are used to calculate predictive properties
#such as Vapor Pressure, Liquid Density, etc. 

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

#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
#This function retrieves experimental data from a large data set to be stored
#for later use. Operates in the same way GetChemicals does.

def experimentalvalues(chemicals):
    SteppChemicalNames = []
    index_list = []

    stepp_data = pd.read_excel('steppdata.xlsx', 'Main Data')
    ChemicalNames = stepp_data['PREFERRED_NAME']

    for i in ChemicalNames:
        SteppChemicalNames.append(i)
        
    for i in chemicals:
        for j in SteppChemicalNames:
            if i==j:
                index_list.append(SteppChemicalNames.index(i))
        
    experimentalvapordata = stepp_data.loc[index_list, ['VAPOR_PRESSURE_MMHG_TEST_PRED']] 

    vapordataframe = pd.DataFrame(experimentalvapordata)

    experimentalvaporvalues = vapordataframe.values

    return experimentalvaporvalues
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
#This function calculates the Liquid Density by using a formula from 
#Physical and Thermodynamic Properties of Pure Chemicals: Data Compilation
#The function takes a list obtained from GetChemicals function

def liquiddensity(lst, t=5):
        
    for i in lst:
        calculated_ld = i/(i*(1-(t/i))**i)
    
    return calculated_ld
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
#This function calculates the vapor pressure by using a formula from 
#Physical and Thermodynamic Properties of Pure Chemicals: Data Compilation
#The function takes a list obtained from GetChemicals function

def vaporpressure(lst, t=5):
    
    for i in lst:
        calculated_vp = np.exp(i + (i/t) + i*np.log(t) + i*t**i)
    
    return calculated_vp
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
#This function creates a database to store all of the theoretical and
#experimental data for n number of chemicals

def properties(chemicals, t=5):
    
    #Empty list to store calculations
    #List of Chemicals
    allproperties = []
    theoretical_vaporpressures = []
    theoretical_liquiddensities = []
    experimental_vaporpressures = []
    experimental_liquiddensities = []         
    chemicallist = chemicals        
    
    #Store variables to get calculations
    chemproperties = GetChemical(chemicals) 
    #Returns a list of chemical properties for each chemical ([[Chemial1_Properties],...,[ChemicalN_Properties]])
    
    #Theoretical
    #Add calculations to list
    for i in chemproperties:
        theoretical_vaporpressures.append(vaporpressure(i))
        theoretical_liquiddensities.append(liquiddensity(i))
    #For each chemical in list, add each vapor pressure to a list
        
    allproperties.append(vaporpressures)
    allproperties.append(experimentalvalues(chemicallist))
    allproperties.append(liquiddensities)
    
    
    #Convert to datagrame
    results = pd.DataFrame(data=allproperties, index=['Predicted Vapor Pressure', 'Experimental Vapor Pressure', 'Predicted Liquid Density'], columns=chemicallist)
                                                                                 
    return results
#-----------------------------------------------------------------------------#



#-----------------------------------------------------------------------------#
#This function saves the dataframe from "properties"

def save(results):
    file = properties(results)
    file.to_csv('Chemical_Properties.csv')    
#-----------------------------------------------------------------------------#

test = properties(['Diethylstilbestrol', 'DDT', 'Ethylenediaminetetraacetic acid'])
print(test)
save(test)