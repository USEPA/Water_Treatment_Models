import csv
from matplotlib import pyplot
import numpy as np
from os.path import isfile

def CSV_Read(Filename):
    with open(Filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count==0:
                Header=row
                lstf=[[] for i in range(len(Header))]
            else:
                index = 0
                for item in row:
                    lstf[index].append(item)
                    index+=1
                line_count += 1
            return [line_count-1, Header, lstf]
                