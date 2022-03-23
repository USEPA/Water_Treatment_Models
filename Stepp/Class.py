from FunctionLibrary import liquiddensity, vaporpressure, GetChemical
import pandas as np
import pandas as pd


class STEPP():
    
    def __init__(self, chemical, t = 15, p = 101325):
        
        TestedChemical = GetChemical(chemical)
        
        #self.function = property of self
        self.calculated_vaporpressure = VapePressure(TestedChemical)
        self.calculated_liquiddensity = LiquidDensity(TestedChemical)
        
        self.properties = pd.DataFrame(data=[], colummns=[])
        

water = STEPP("Water")
print(water.properties)