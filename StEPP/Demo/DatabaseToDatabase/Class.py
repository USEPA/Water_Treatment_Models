from FunctionLibrary import LiquidDensity, VapePressure, GetChemical
import pandas as np
import pandas as pd


class STEPP():
    
    def __init__(self, chemical, t=15, p = 101325):
     
        
        #self.function = property of self
        self.VaporPressure = VapePressure(data)
        self.LiquidDensity = LiquidDensity(data)
        
        
data = STEPP(GetChemical("Water"))
        
        
        
        
