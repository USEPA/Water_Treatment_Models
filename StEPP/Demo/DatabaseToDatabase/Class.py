from FunctionLibrary import LiquidDensity, VapePressure
import pandas as np
import pandas as pd

#WIP

class STEPP():
    
    def __init__(self, chemical, t=15, p = 101325):
     
        #Read in file
        ChemicalData = pd.read_csv('Chemicals.csv')
     
        #Get constants
        #We want to index multiple constants of a chemical
        #We want to index multiple columns of a row
        
        ChemicalVariableData = Chemicaldata.loc[1, ['A','B','C','D','E']]
        
        #Create Dataframe
        #Put constants in dataframe_constants
        
        ChemicalConstants = pd.DataFrame(ChemicalVariableData)
        
        #Assign data for function
        FormulaConstants = ChemicalVariableData.values
        
        #self.function = property of self
        self.vaporpressure = VaporPressure(FormulaConstants)
            
        
        
        
        
        
