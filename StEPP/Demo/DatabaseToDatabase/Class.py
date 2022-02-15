import pandas as pd
import numpy as np

class Chemical(object):
    
    def __init__(self, database):
            
            self.database = database
            data = pd.read_csv('Chemicals.csv')
            data.set_index("Property", inplace = True)
            
    def GetLDVariables(self, A, B, C, D, E, T):
        
        rowA = data.loc[A]                                  #Get A row
        rowB = data.loc[B]                                  #Get B row
        rowC = data.loc[C]                                  #Get C row
        rowD = data.loc[D]                                  #Get D row
        rowT = data.loc[T]                                  #Get T row

        Data = []                                         #Create an Empty list to store data

        Data.append(rowA)                                 #Add row A
        Data.append(rowB)                                 #Add row B
        Data.append(rowC)                                 #Add row C
        Data.append(rowD)                                 #Add row D
        Data.append(rowT)                                 #Add row T

        Database = pd.DataFrame(Data)                     #Create a data frame
        Database.to_csv('ChemicalOutput.csv', index = False) #Create CSV file out of Data frame


    Vars = GetLDVariables("A", "B", "C", "D", "E", "T")   #Get wanted variables
    ab = pd.read_csv('ChemicalOutput.csv')                #Read new CSV file
        

    def LiquidDensity(self, Chemical):
       
        ChemicalProperties = ab[Chemical]
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
    
    def GetVPVariables(A, B, C, D, E, T):
        
        try:
            df = pd.read_csv('Chemicals.csv')                 #ReadCSVDatabase

            df.set_index("Property", inplace=True)

            Chemical = df.columns                             #Get Header
        
            rowA = df.loc[A]                                  #Get A row
            rowB = df.loc[B]                                  #Get B row
            rowC = df.loc[C]                                  #Get C row
            rowD = df.loc[D]                                  #Get D row
            rowE = df.loc[E]                                  #Get E row
            rowT = df.loc[T]                                  #Get T row

            Data = []                                         #Create an Empty list to store data

            Data.append(rowA)                                 #Add row A
            Data.append(rowB)                                 #Add row B
            Data.append(rowC)                                 #Add row C
            Data.append(rowD)                                 #Add row D
            Data.append(rowE)                                 #Add row E
            Data.append(rowT)                                 #Add row T

            Database = pd.DataFrame(Data)                     #Create a data frame
            Database.to_csv('ChemicalOutput.csv', index = False) #Create CSV file out of Data frame

        except KeyError:
            print('Data not found in database')
            
    Vars = GetLDVariables("A", "B", "C", "D", "E", "T")

    def VapePressure(Chemical):
        
        ab = pd.read_csv('ChemicalOutput.csv')            #Read new CSV file
        
        ChemicalProperties = ab[Chemical]
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
    
    def GetVPVariables(Mass):
        
        try:
            df = pd.read_csv('Chemicals.csv')                 #ReadCSVDatabase

            df.set_index("Property", inplace=True)

            Chemical = df.columns                             #Get Header
        
            rowMass = df.loc[Mass]                                  #Get Mass row
          
            Data = []                                         #Create an Empty list to store data

            Data.append(rowMass)                                 #Add row Mass
           
            Database = pd.DataFrame(Data)                     #Create a data frame
            Database.to_csv('ChemicalOutput.csv', index = False) #Create CSV file out of Data frame

        except KeyError:
            print('Data not found in database')
            
    Vars = GetLDVariables("A", "B", "C", "D", "E", "T")

    def MolarVolume(Chemical):
        
        ChemicalProperties = df[Chemical]
        #Column for Specific Chemical
        
        ChemicalMass = ChemicalProperties[0]
        LD = LiquidDensity.LiquidDensity(Chemical)
        
        Volume = ChemicalMass/LiquidDensity
        #Formula from 
        #Physical and Thermodynamic Properties of Pure Chemicals: Data Compilation

        return Volume

    print(MolarVolume("Water"))
    
    def LiquidDiffusivity(Chemical):
        
        mu = 0.00131 * 1000
        #Viscosity of Water in g/cm/s to cp
        
        Diffusivity = 13.26/((mu)**1.4*MolarVolume(Chemical)*100**2)
        #Prediction of Diffusion Coefficients for Non-electrolytes in Dilute 
        #Aqueous Solutions
        
        return Diffusivity

        