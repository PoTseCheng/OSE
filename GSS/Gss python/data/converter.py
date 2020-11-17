import scipy.io
import numpy as np
import os



def converter():
    '''
    convert .mat to csv files
    '''
    dirname=os.getcwd()
    dirname=dirname + "/data/"
    data = scipy.io.loadmat(dirname+"epsi_test.mat")
    for i in data:
        if '__' not in i and 'readme' not in i: 
            np.savetxt("epsi_test.csv",data[i],delimiter=',')
    data = scipy.io.loadmat(dirname+"epsi10000.mat")
    for i in data:
        if '__' not in i and 'readme' not in i: 
            np.savetxt("epsi10000.csv",data[i],delimiter=',')


            
    


