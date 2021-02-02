from numpy import sqrt
from numpy.random import rand
from scipy.special import erfinv
import numpy as np
import math
from scipy import linalg

#Borrow function from others to ensure same results of function Randn
def randn2(*args,**kwargs):
    '''
    Calls rand and applies inverse transform sampling to the output. Borrowed from Jonas Rauber.
    '''
    uniform = rand(*args, **kwargs)
    return sqrt(2) * erfinv(2 * uniform - 1)
    # Copyright (c) 2015 Jonas Rauber

def Productivity(T, N, a_init, sigma, rho):
    '''
    This function simulates draws of random series of the productivity shocks and the corresponding series of the productivity levels.
    ----------
    Arguments:
    T(int): Simulation length, needs to be at least 1
    N(int): Number of countries, notice N needs to be at least 1.
    a_int(2D-numpy): The initial condition given for the productivity levels of N countries.
    sigma(float): Parameters of the model in Kenneth L. Judd et.al (2011).
    rho(float): Parameters of the model in Kenneth L. Judd et.al (2011).

    ----------
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
    ----------
    Arguments:
    z(2D numpy array): Data points on which the polynomial basis functions must be constructed.
    D(int): 1 to 5 the degree of the polynomial whose basis functions must be constructed.

    ----------
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

def GH_Quadrature(Qn, N, vcv):
    '''
    This function constructs integration nodes and weights under Gauss-Hermite quadrature (product) integration rule with Qn<=10 nodes in each of N dimensions.
    ----------
    Arguments:
    Qn(int): The number of nodes in each dimension, notice that Qn must be between 1 to 10.
    N(int): The number of countries.
    vcv(2D numpy array): The predefined covariance matrix.

    ----------
    Outputs:
    n_nodes(int):Total number of integration nodes.
    epsi_nodes(2D numpy array):Integration nodes.
    weight_nodes(2D numpy array):Integration weights.
    
    '''
    #For the following codes#######################
    #Qn :Number of nodes in each dimension; Qn <=10
    #eps: Set of integration nodes
    #weight. Set of integration weights      
    if Qn == 1:
        eps = [0]
        weight = [math.sqrt(math.pi)]
    elif Qn == 2:
        eps = [0.7071067811865475, -0.7071067811865475]
        weight = [0.8862269254527580,  0.8862269254527580]
    elif Qn == 3:
        eps = [1.224744871391589, 0, -1.224744871391589]
        weight = [0.2954089751509193, 1.181635900603677, 0.2954089751509193]
    elif Qn == 4:
        eps = [1.650680123885785, 0.5246476232752903, -0.5246476232752903, -1.650680123885785]
        weight = [0.08131283544724518, 0.8049140900055128, 0.8049140900055128, 0.08131283544724518]
    elif Qn == 5:
        eps = [2.020182870456086, 0.9585724646138185, 0, -0.9585724646138185, -2.020182870456086]
        weight = [0.01995324205904591, 0.3936193231522412, 0.9453087204829419, 0.3936193231522412, 0.01995324205904591]
    elif Qn == 6:
        eps = [2.350604973674492, 1.335849074013697, 0.4360774119276165, -0.4360774119276165, -1.335849074013697, -2.350604973674492]
        weight = [0.004530009905508846, 0.1570673203228566, 0.7246295952243925, 0.7246295952243925, 0.1570673203228566, 0.004530009905508846]
    elif Qn == 7:
        eps = [2.651961356835233, 1.673551628767471, 0.8162878828589647, 0, -0.8162878828589647, -1.673551628767471, -2.651961356835233]
        weight = [0.0009717812450995192, 0.05451558281912703, 0.4256072526101278, 0.8102646175568073, 0.4256072526101278, 0.05451558281912703, 0.0009717812450995192]
    elif Qn == 8:
        eps = [2.930637420257244, 1.981656756695843, 1.157193712446780, 0.3811869902073221, -0.3811869902073221, -1.157193712446780, -1.981656756695843, -2.930637420257244]
        weight = [0.0001996040722113676, 0.01707798300741348, 0.2078023258148919, 0.6611470125582413, 0.6611470125582413, 0.2078023258148919, 0.01707798300741348, 0.0001996040722113676]
    elif Qn == 9:
        eps = [3.190993201781528, 2.266580584531843, 1.468553289216668, 0.7235510187528376, 0, -0.7235510187528376, -1.468553289216668, -2.266580584531843, -3.190993201781528]
        weight = [0.00003960697726326438, 0.004943624275536947, 0.08847452739437657, 0.4326515590025558, 0.7202352156060510, 0.4326515590025558, 0.08847452739437657, 0.004943624275536947, 0.00003960697726326438]
    else:
        Qn = 10
        eps = [3.436159118837738, 2.532731674232790, 1.756683649299882, 1.036610829789514, 0.3429013272237046, -0.3429013272237046, -1.036610829789514, -1.756683649299882, -2.532731674232790, -3.436159118837738]
        weight = [7.640432855232621e-06, 0.001343645746781233, 0.03387439445548106, 0.2401386110823147, 0.6108626337353258, 0.6108626337353258, 0.2401386110823147, 0.03387439445548106, 0.001343645746781233, 7.640432855232621e-06]
    #N-dimensional integration nodes and weights for N uncorrelated normally distributed random variables with zero mean and unit variance
    
    #Total number of integration nodes
    n_nodes = Qn**N
    
    #A better approch for 2D array construction
    z1 = np.ones([n_nodes,N]).astype(float)
    w1i = np.ones([n_nodes,N]).astype(float)
    w1 = np.ones([n_nodes,1]).astype(float)
    for i in range(N):
        z1[:,i]=np.tile(np.repeat(eps,10**(i)), 10**(N-i-1))
        w1i[:,i]=np.tile(np.repeat(weight,10**(i)), 10**(N-i-1))
    for i in range(N):
        w1[:,0]*=w1i[:,i]
    
    #Integration nodes preparations
    z = math.sqrt(2)*z1
    #Integration weights are the same for the cases of correlated and uncorrelated random variables 
    weight_nodes = w1/(math.sqrt(math.pi)**N)
    
    sqrt_vcv = linalg.cholesky(vcv)
    #Final integration nodes
    epsi_nodes = z@sqrt_vcv

    return n_nodes, epsi_nodes, weight_nodes

def Monomials_1(N, vcv):
    '''
    This function constructs integration nodes and weights under N-dimensional monomial (non-product) integration rule with 2N nodes.
    ----------
    Arguments:
    N(int): Number of countries.
    vcv(2D numpy array): The predefined covariance matrix.
    
    ----------
    Outputs:
    n_nodes(int): Total number of integration nodes.
    epsi_nodes(2D numpy array): Integration nodes.
    weight_nodes(2D numpy array): Integration weights.
    
    '''

    #Total number of integration nodes
    n_nodes = 2*N

    #Step 1. N-dimensional integration nodes for N uncorrelated random variables with zero mean and unit variance
    ####################################################################################################

    #construct container for values
    z1 = np.zeros((n_nodes, N))

    #In each node, random variable i takes value either 1 or -1, and all other variables take value 0
    for i in range(1, N+1):
        z1[2*(i-1):2*i,i-1]= np.array([1,-1])

    #Step 2. N-dimensional integration nodes and weights for N correlated random variables with zero mean and variance-covaraince matrix vcv
    #####################################################################################################

    #Preparations
    sqrt_vcv = linalg.cholesky(vcv)
    R = math.sqrt(N)*sqrt_vcv

    #Integration nodes
    epsi_nodes = z1@R
    
    #Integration weights
    weight_nodes = np.ones((n_nodes,1))/n_nodes

    return n_nodes, epsi_nodes, weight_nodes


def Monomials_2(N, vcv):
    '''
    This function constructs integration nodes and weights under N-dimensional monomial (non-product) integration rule with 2N^2+1 nodes.
    ----------
    Arguments:
    N(int): Number of countries.
    vcv(2D numpy array): The predefined covariance matrix.

    ----------
    Outputs:
    n_nodes(int): Total number of integration nodes.
    epsi_nodes(2D numpy array): Integration nodes.
    weight_nodes(2D numpy array): Integration weights.
    '''
    #Total number of integration nodes
    n_nodes = 2*N**2+1

    #Step 1: N-dimensional integration nodes for N uncorrelated random variables with zero mean and unit variance
    ####################################################################################################

    #Point origin(0 dimension)
    z0 = np.zeros((1,N))

    #Deviations in one dimension
    #In each node, random variable i takes value either 1 or -1, and all other variables take value 0
    z1 = np.zeros((2*N,N))
    for i in range(1, N+1):
        z1[2*(i-1):2*i,i-1]= np.array([1,-1])
    
    #Deviations in 2nd dimension
    #In each node, a pair of random variables (p,q) takes either values (1,1) or (1,-1) or (-1,1) or (-1,-1), and all other variables take value 0

    z2 = np.zeros((2*N*(N-1),N))
    i = 0
    for p in range(1, N):
        for q in range(p+1, N+1):
            i +=1
            z2[4*(i-1):4*i,p-1] = np.array([1, -1, 1, -1])
            z2[4*(i-1):4*i,q-1] = np.array([1, 1, -1, -1])

    # Step2: N-dimensional integration nodes and weights for N correlated random variables with zero mean and variance-covaraince matrix vcv
    ########################################################################################################################################

    #Preparations
    sqrt_vcv = linalg.cholesky(vcv)
    R = math.sqrt(N+2)*sqrt_vcv
    S = math.sqrt((N+2)/2)*sqrt_vcv

    #Integration nodes
    epsi_nodes = np.vstack((z0, z1@R, z2@S))

    #Integration weights
    #See condition in (B.8) in the Supplement of Judd et al. (2011) 
    weight_nodes = np.vstack((
                        2/(N+2)*np.ones((z0.shape[0], 1)),
                        (4-N)/2/(N+2)**2*np.ones((z1.shape[0], 1)),
                        1/(N+2)**2*np.ones((z2.shape[0], 1))
                    ))
    
    return n_nodes, epsi_nodes, weight_nodes

