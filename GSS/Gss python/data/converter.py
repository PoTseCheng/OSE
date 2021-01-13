import scipy.io
import numpy as np
import os


def converter_new(x):
    '''
    input:
          x(.mat): mat file in data folder
    output:
          a csv file at the main folder
    '''
    dirname = os.getcwd()
    dirname = dirname + "/data/"
    data = scipy.io.loadmat(dirname+x)
    
    for i in data:
        if '__' not in i and 'readme' not in i: 
            np.savetxt( x + ".csv", data[i], delimiter=',')

    return