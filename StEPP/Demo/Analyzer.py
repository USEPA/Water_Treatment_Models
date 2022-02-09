from CSV_Reader import CSV_Read

class Analyzer(object):
    def __init__(self, database):
        # Load Database
        self.database = database
        csv = CSV_Read(database)
        self.number_of_datasets = csv[0]
        self.datasets = csv[2]
        self.chemical_index = csv[1].index("Chemicals")


    def Reciever(self, dataset_index):
        # Get information about Chemical
        value_col = self.datasets[self.chemical_index][dataset_index]

        # Get index
        data_chemical_index = header.index(value_col)
        
