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
    This function simulates draws of random series of the productivity shocks and the corresponding series of the productivity levels.
    Arguments:
    T(int): Simulation length, needs to be at least 1
    N(int): Number of countries, notice N needs to be at least 1.
    a_int(2D-numpy): The initial condition given for the productivity levels of N countries.
    sigma(float): Parameters of the model in Kenneth L. Judd et.al (2011).
    rho(float): Parameters of the model in Kenneth L. Judd et.al (2011).

    Output:
    a(2D numpy array): Time series of the productivity levels of N countries.
    
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
        a = np.concatenate((a,a[i]**rho*np.exp(epsi[i+1]).reshape(1,N)), axis=0) 
    
    return a

def Ord_Polynomial_N(z, D):
    '''
    Summary
    Arguments:
    z(2D numpy array): Data points on which the polynomial basis functions must be constructed.
    D(int): 1 to 5 the degree of the polynomial whose basis functions must be constructed.

    Output:
    basis_fs(2D numpy array): Matrix with complete degrees of polynominal.
    
    '''
    n_rows, dimen = np.shape(z)
    # The matrix of the basis functions of the first-degree polynomial (Default)
    basis_fs = np.hstack((np.ones((n_rows, 1)),z))
    #  Index number of a polynomial basis function
    
    #2nd-degree polynomial:
    if D == 2:
        for i in range(n_rows):
            tempt = (z[i][:,np.newaxis]@z[i][np.newaxis,:])[np.tril_indices(dimen)]
            basis_fs= np.hstack((basis_fs[i], tempt))[np.newaxis,]

    #3rd-degree polynomial:
    elif D == 3:
        for i in range(1, dimen+1):
            for j in range(i, dimen+1):
                basis_fs = np.hstack((basis_fs, z[:,i-1]*z[:,j-1][np.newaxis,]))
                for k in range(j, dimen+1):
                    basis_fs = np.hstack((basis_fs, z[:,i-1]*z[:,j-1]*z[:,k-1][np.newaxis,]))
    #4th-degree polynominal:
    elif D == 4:
        for i in range(1, dimen+1):
            for j in range(i, dimen+1):
                basis_fs = np.hstack((basis_fs, z[:,i-1]*z[:,j-1][np.newaxis,]))
                for k in range(j, dimen+1):
                    basis_fs = np.hstack((basis_fs, z[:,i-1]*z[:,j-1]*z[:,k-1][np.newaxis,]))
                    for l in range(k, dimen+1):
                        basis_fs = np.hstack((basis_fs, z[:,i-1]*z[:,j-1]*z[:,k-1]*z[:,l-1][np.newaxis,]))
    #5th-degree polynominal:
    elif D == 5:
        for i in range(1, dimen+1):
            for j in range(i, dimen+1):
                basis_fs = np.hstack((basis_fs, z[:,i-1]*z[:,j-1][np.newaxis,]))
                for k in range(j, dimen+1):
                    basis_fs = np.hstack((basis_fs, z[:,i-1]*z[:,j-1]*z[:,k-1][np.newaxis,]))
                    for l in range(k, dimen+1):
                        basis_fs = np.hstack((basis_fs, z[:,i-1]*z[:,j-1]*z[:,k-1]*z[:,l-1][np.newaxis,]))
                        for m in range(l, dimen+1):
                            basis_fs = np.hstack((basis_fs, z[:,i-1]*z[:,j-1]*z[:,k-1]*z[:,l-1]*z[:,m-1][np.newaxis,]))

    return basis_fs


