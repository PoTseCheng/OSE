import scipy.io
import numpy as np
import os



def converter():
    '''
    convert .mat to csv files
    '''
    dirname=os.getcwd()
    dirname=dirname + "/data/"
    data = scipy.io.loadmat(dirname+"DataAgriculture.mat")
    for i in data:
        if '__' not in i and 'readme' not in i: 
            np.savetxt("DataAgriculture.csv",data[i],delimiter=',')


            
    


