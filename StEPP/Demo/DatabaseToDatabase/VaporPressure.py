import pandas as pd
import numpy as np

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