import pandas as pd
import numpy as np
import LiquidDensity

df = pd.read_csv('Chemicals.csv', header = 0)
#Load CSV File

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
