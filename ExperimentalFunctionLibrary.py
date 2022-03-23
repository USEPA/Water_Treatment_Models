import numpy as np
import pandas as pd


#-----------------------------------------------------------------------------#
#This function retrieves experimental data from a large data set to be stored
#for later use. Operates in the same way GetChemicals does.

def experimentalvalues(chemicals):
    
    #Define Empty Lists for Later Use
    SteppChemicalNames = []             #Stores all available chemical names
    index_list = []                     #Row numbers we want
    experimentaldata = []               #Stores all experimental data we want
    molarvolume = []                    #Molar Volume is a calculation
    chemicallist = chemicals

    stepp_data = pd.read_excel('SteppChemicals.xlsx', 'Sheet1')
    ChemicalNames = stepp_data['PREFERRED_NAME']


    #List of all chemical names
    for i in ChemicalNames:
        SteppChemicalNames.append(i)
        
    #For element in chemicals, if element is in ChemicalList (lst), get index
    #of that element in ChemicalList    
    for i in chemicals:
        for j in SteppChemicalNames:
            if i==j:
                index_list.append(SteppChemicalNames.index(i))
        
    #Retrieve wanted experimental data    
    experimentalvapordata = stepp_data.loc[index_list, ['VAPOR_PRESSURE_MMHG_TEST_PRED']] 
    experimentalliquiddensity = stepp_data.loc[index_list, ['DENSITY_G/CM^3_TEST_PRED']]
    experimentalmolecularweight = stepp_data.loc[index_list, ['AVERAGE_MASS']]
    experimentalboilingpoint = stepp_data.loc[index_list, ['BOILING_POINT_DEGC_TEST_PRED']]
    experimentalsolubility = stepp_data.loc[index_list, ['WATER_SOLUBILITY_MOL/L_TEST_PRED']]
    experimentalkoa = stepp_data.loc[index_list, ['OCTANOL_AIR_PARTITION_COEFF_LOGKOA_OPERA_PRED']]
    
    
    # #We want to put the vapor and liquid values into an 'experimental data frame'
    evd=experimentalvapordata.values
    eld=experimentalliquiddensity.values
    emw=experimentalmolecularweight.values
    ebp=experimentalboilingpoint.values
    es=experimentalsolubility.values
    ekoa=experimentalkoa.values
    
    #With the way the data is stored by default I flattened to make my life easier
    flatevd=list(np.concatenate(evd).flat)
    flateld=list(np.concatenate(eld).flat)
    flatemw=list(np.concatenate(emw).flat)
    flatebp=list(np.concatenate(ebp).flat)
    flates=list(np.concatenate(es).flat)
    flatekoa=list(np.concatenate(ekoa).flat)
    
    #Caluclation for molar volume
    for mass,density in zip(flatemw,flateld):
            if type(mass) and type(density) == float:
                molarvolume.append(mass/density)
            else:
                molarvolume.append("-")
    
    #Store all the data for a dataframe
    experimentaldata.append(flatevd)
    experimentaldata.append(flateld)
    experimentaldata.append(flatemw)
    experimentaldata.append(flatebp)
    experimentaldata.append(flates)
    experimentaldata.append(flatekoa)
    experimentaldata.append(molarvolume)

    return experimentaldata

print(experimentalvalues(['Diethylstilbestrol', 'DDT', 'Ethylenediaminetetraacetic acid', 'Hydrogen']))     
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

    #Returns a list of chemical properties for each chemical ([[Chemial1_Properties],...,[ChemicalN_Properties]])    
        
    for i in experimentalvalues(chemicallist):
        allproperties.append(i)
    
    
    #Convert to datagrame
    results = pd.DataFrame(data=allproperties, index=['Measured Vapor Pressure', 'Measured Liquid Density', 'Measured Molecular Weight', 'Measured Boiling Point', 'Measured Water Solubility', 'Measured Octonal Air Partition Coefficient', 'Measured Volume'], columns=chemicallist)
                                                                                 
    return results
    
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
        
    allproperties.append(theoretical_vaporpressures)
    allproperties.append(theoretical_liquiddensities)
    allproperties.append(experimentalvalues(chemicallist))
    
    #Convert to datagrame
    results = pd.DataFrame(data=allproperties, index=['Predicted Vapor Pressure', 'Predicted Liquid Density',], columns=chemicallist)
                                                                                 
    return results
#-----------------------------------------------------------------------------#




#-----------------------------------------------------------------------------#
#This function saves the dataframe from "properties"

def save(props):
    file = props
    file.to_csv('Chemical_Properties3.csv')    
#-----------------------------------------------------------------------------#




testchems=['Diethylstilbestrol', 'DDT', 'Ethylenediaminetetraacetic acid', 'Hydrogen']
test=properties(testchems)
print(test)
save(test)