from numpy import sqrt
from numpy.random import rand
from scipy.special import erfinv
import numpy as np


#Borrow function from others to ensure same results of function Randn
def randn2(*args,**kwargs):
    '''
    Calls rand and applies inverse transform sampling to the output.
    '''
    uniform = rand(*args, **kwargs)
    return sqrt(2) * erfinv(2 * uniform - 1)
    # Copyright (c) 2015 Jonas Rauber
    # License: The MIT License (MIT)
    # See: https://github.com/jonasrauber/randn-matlab-python

def Productivity(T, N, a_init, sigma, rho):
    '''
    Summary
    Arguments:
    T(int):
    N(int):
    a_int(2D-numpy):
    sigma(float):
    rho(float):

    Output:
    a(2D numpy array):
    
    '''
    np.random.seed(123)
    #random draw of common-for-all-countries productivity shocks for T periods
    EPSI = randn2(T,1)
    #random draw of country-specific productivity shocks for T periods and N countries
    epsi = randn2(T,N)
    #Compute the error terms in the process for productivity level using condition (4) in JMM (2011)
    epsi = sigma*(epsi+ EPSI@np.ones((1,N)))
    #Initial condition for the productivity levels
    a = a_init

    for i in range(T-1):
        a = np.concatenate((a,a[i]**rho*np.exp(epsi[i+1]).reshape(1,N)), axis =0) 
    
    return a