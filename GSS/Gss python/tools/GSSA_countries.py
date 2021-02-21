from numpy import sqrt
from numpy.random import rand
from scipy.special import erfinv
import numpy as np
import math
from scipy import linalg
from scipy.optimize import linprog
import scipy.io
import time
import pandas as pd
import os

def randn2(*args,**kwargs):
    '''
    The randn function that matlab uses.
    '''
    uniform = rand(*args, **kwargs)
    return sqrt(2) * erfinv(2 * uniform - 1)
 

def Productivity(T, N, a_init, sigma, rho):
    '''
    This function simulates draws of random series of the productivity shocks 
    and the corresponding series of the productivity levels.
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
        for i in range(1, dimen+1):
            for j in range(i, dimen+1):
                basis_fs = np.hstack((basis_fs, (z[:,i-1]*z[:,j-1]).reshape(n_rows,1)))
        #for i in range(n_rows):
            #basis_fs= np.hstack((basis_fs[i],
            #(z[i][:,np.newaxis]@z[i][np.newaxis,:])[np.tril_indices(dimen)]
                                #))[np.newaxis,]

    #3rd-degree polynomial:
    elif D == 3:
        for i in range(1, dimen+1):
            for j in range(i, dimen+1):
                basis_fs = np.hstack((basis_fs, (z[:,i-1]*z[:,j-1]).reshape(n_rows,1)))
                for k in range(j, dimen+1):
                    basis_fs = np.hstack((basis_fs, (z[:,i-1]*z[:,j-1]*z[:,k-1]).reshape(n_rows,1)))
    #4th-degree polynominal:
    elif D == 4:
        for i in range(1, dimen+1):
            for j in range(i, dimen+1):
                basis_fs = np.hstack((basis_fs, (z[:,i-1]*z[:,j-1]).reshape(n_rows,1)))
                for k in range(j, dimen+1):
                    basis_fs = np.hstack((basis_fs, (z[:,i-1]*z[:,j-1]*z[:,k-1]).reshape(n_rows,1)))
                    for l in range(k, dimen+1):
                        basis_fs = np.hstack((basis_fs, (z[:,i-1]*z[:,j-1]*z[:,k-1]*z[:,l-1]).reshape(n_rows,1)))
    #5th-degree polynominal:
    elif D == 5:
        for i in range(1, dimen+1):
            for j in range(i, dimen+1):
                basis_fs = np.hstack((basis_fs, (z[:,i-1]*z[:,j-1]).reshape(n_rows,1)))
                for k in range(j, dimen+1):
                    basis_fs = np.hstack((basis_fs, (z[:,i-1]*z[:,j-1]*z[:,k-1]).reshape(n_rows,1)))
                    for l in range(k, dimen+1):
                        basis_fs = np.hstack((basis_fs, (z[:,i-1]*z[:,j-1]*z[:,k-1]*z[:,l-1]).reshape(n_rows,1)))
                        for m in range(l, dimen+1):
                            basis_fs = np.hstack((basis_fs, (z[:,i-1]*z[:,j-1]*z[:,k-1]*z[:,l-1]*z[:,m-1]).reshape(n_rows,1)))

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
        eps = np.zeros([1,1])
        weight = math.sqrt(math.pi)
    elif Qn == 2:
        eps = np.array([0.7071067811865475, -0.7071067811865475])
        weight = np.array([0.8862269254527580,  0.8862269254527580])
    elif Qn == 3:
        eps = np.array([1.224744871391589, 0., -1.224744871391589])
        weight = np.array([0.2954089751509193, 1.181635900603677, 0.2954089751509193])
    elif Qn == 4:
        eps = np.array([1.650680123885785, 0.5246476232752903, -0.5246476232752903, -1.650680123885785])
        weight = np.array([0.08131283544724518, 0.8049140900055128, 0.8049140900055128, 0.08131283544724518])
    elif Qn == 5:
        eps = np.array([2.020182870456086, 0.9585724646138185, 0.,-0.9585724646138185,-2.020182870456086])
        weight = np.array([0.01995324205904591,0.3936193231522412,0.9453087204829419,0.3936193231522412,0.01995324205904591])
    elif Qn == 6:
        eps = np.array([2.350604973674492,1.335849074013697,0.4360774119276165,-0.4360774119276165,-1.335849074013697,-2.350604973674492])
        weight = np.array([0.004530009905508846,0.1570673203228566,0.7246295952243925,0.7246295952243925,0.1570673203228566,0.004530009905508846])
    elif Qn == 7:
        eps = np.array([2.651961356835233,1.673551628767471,0.8162878828589647,0.,-0.8162878828589647,-1.673551628767471,-2.651961356835233])
        weight = np.array([0.0009717812450995192, 0.05451558281912703,0.4256072526101278,0.8102646175568073,0.4256072526101278,0.05451558281912703,0.0009717812450995192])
    elif Qn == 8:
        eps = np.array([2.930637420257244,1.981656756695843,1.157193712446780,0.3811869902073221,-0.3811869902073221,-1.157193712446780,-1.981656756695843,-2.930637420257244])
        weight = np.array([0.0001996040722113676,0.01707798300741348,0.2078023258148919,0.6611470125582413,0.6611470125582413,0.2078023258148919,0.01707798300741348,0.0001996040722113676])
    elif Qn == 9:
        eps = np.array([3.190993201781528,2.266580584531843,1.468553289216668,0.7235510187528376,0,-0.7235510187528376,-1.468553289216668,-2.266580584531843,-3.190993201781528])
        weight = np.array([0.00003960697726326438,0.004943624275536947,0.08847452739437657,0.4326515590025558,0.7202352156060510,0.4326515590025558,0.08847452739437657,0.004943624275536947,0.00003960697726326438])
    else:
        Qn = 10
        eps = np.array([3.436159118837738,2.532731674232790,1.756683649299882,1.036610829789514,0.3429013272237046,-0.3429013272237046,-1.036610829789514,-1.756683649299882,-2.532731674232790,-3.436159118837738])
        weight = np.array([7.640432855232621e-06,0.001343645746781233,0.03387439445548106,0.2401386110823147,0.6108626337353258,0.6108626337353258,0.2401386110823147,0.03387439445548106,0.001343645746781233,7.640432855232621e-06])
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

def Num_Stab_Approx(X, Y, RM, penalty, normalize):
    '''
    This function implements the approximation methods mentioned in Judd et al. (2011)
    --------
    Arguments:
        X(2D numpy array): Matrix of dependent variables in a regression.
        Y(2D numpy array):Matrix of independent variables.
        RM(int): Regression (approximation) method, from 1 to 6. RLAD-DP & LAD-DP are not included, see notebook for reasoning.
        penalty(int): Regularisation parameter for a regularisation methods.
        normalize(int): Optional parameter to normalise the data or not. 1 or 0.

    ---------
    Outputs:
        B(2D numpy array): Matrix of the regression coefficients.

    '''
    
    # Step 1: Compute the dimensionality of the data
    ################################################

    T, n = X.shape
    N  = Y.shape[1]
    
    # Step 2: Normalize the data(or not)
    ####################################

    if ((normalize==1) or (RM>=5)):
        X1 = (X[:,1:n]-np.ones((T,1))@X[:,1:n].mean(axis=0)[np.newaxis,])/(np.ones((T,1))@np.std(X[:,1:n], axis=0, ddof=1)[np.newaxis,])
        Y1 = (Y - np.ones((T,1))@Y.mean(axis=0)[np.newaxis,])/(np.ones((T,1))@np.std(Y, axis=0, ddof= 1)[np.newaxis,])
        n1 = n-1
    
    else:#leave the values unchange if not normalised
        X1 = X
        Y1 = Y
        n1 = n

    # Step 3: Regression methods
    ############################

    #OLS
    if RM == 1:
        B = linalg.inv(X1.conj().T@X1)@X1.conj().T@Y1

    #LS-SVD
    elif RM == 2:
        U, S, Vh = linalg.svd(X1 ,full_matrices=False)
        V = Vh.T
        S_inv = np.diag(1/S)
        B = V@S_inv@U.conj().T@Y1

    #LAD-PP
    elif RM == 3:
        BND = [(-100, 100)]*n1 + [(0, None)]*2*T
        f = np.vstack((np.zeros((n1,1)), np.ones((2*T,1))))
        Aeq = np.concatenate((X1, np.eye(T), -np.eye(T)), axis=1)
        B =[]
        #solve the equation
        for i in range(N):
            beq = Y1[:,i]
            result = linprog(f, A_eq = Aeq, b_eq = beq, bounds= BND, method="highs-ipm")
            B.append(list(result.x[0:n1]))
        B = np.asarray(B).T

    # RLS-Tikhonov
    elif RM == 4:
        B = linalg.inv(X1.conj().T@X1+T/n1*np.eye(n1)*10**penalty)@X1.conj().T@Y1

    # RLS-TSVD
    elif RM == 5:
        U, S, Vh = linalg.svd(X1, full_matrices=False)
        V = Vh.T
        r = np.count_nonzero(np.divide(np.diag(S).max(), np.diag(S))<= 10**(penalty))
        Sr_inv = np.zeros((n1,n1))
        Sr_inv[0:r, 0:r]= np.diag(np.divide(1., S[:r]))
        B = V@Sr_inv@U.conj().T@Y1

    # RLAD-PP
    elif RM == 6:
        #we can just use the default setting from scipy as the lower and upper will be the same
        f= np.vstack((10**penalty*np.ones((n1*2,1))*T/n1, np.ones((2*T,1))))
        Aeq= np.c_[X1,-X1, np.eye(T), -np.eye(T)]
        B = []
        #solve the equation
        for i in range(N):
            beq = Y1[:,i]
            result = linprog(f, A_eq = Aeq, b_eq = beq, method="highs-ipm")
            B.append(list(result.x[0:n1]-result.x[n1:2*n1]))
        B = np.asarray(B).T
    
    #Step 4: Infer the regression coefficients in the original regression with unnormalised data
    ############################################################################################

    if ((normalize==1) or (RM>=5)):
        B2 = (1/np.std(X[:,1:n], axis=0).conj().T).reshape((n1,1))@np.std(Y, axis=0)[np.newaxis,]*B
        B1 = Y.mean(axis=0)[np.newaxis,] - X[:, 1:n].mean(axis=0)[np.newaxis,]@B2
        B = np.vstack((B1,B2))
    
    return B

def Accuracy_Test_N(k, a, bk, D, IM, alpha, gam, delta, beta, A, tau, rho, vcv, discard):
    '''
    This function is used for evaluating accuracy of solutions to the multi-country model: it computes approximation errors in 
    the optimality conditions on a given set of points in the state space.
    --------
    Arguments:

        k(2D numpy array): Current-period capital.
        a(2D numpy array): Current productivity levels.
        bk(2D numpy array): Coefficients of the capital policy functions of N countries.
        IM(int): Integration method in the original GSSA model.
        alpha(float): Capital share in output.
        gam(float): Utility-function parameter.
        delta(float): Depreciation rate.
        beta(float): Discount factor.
        A(float): The normalizing constant in output.
        tau(float): The welfare weight of country.
        rho(float): Persistence of the log of the productivity level.
        vcv(2D numpy array): Variance-covariance matrix of the countries' productivity shocks.
        discard(int): Data points to discard.
    --------
    Output:
    
        Errors_mean(float): The mean approximation errors.
        Errors_max(float): The maximum approximation errors.
        time_test(float): The time to run the test.
    '''
    start = time.time()
    #1. Get number of points P, which accuracy are evaluated with respect to country N
    P, N = a.shape
    
    #2. Integration method for evaluating accuracy 
    if ((IM>=1) and (IM<=10)):
        n_nodes, epsi_nodes, weight_nodes = GH_Quadrature(IM, N, vcv)
    elif IM == 11:
        n_nodes, epsi_nodes, weight_nodes = Monomials_1(N,vcv)
    elif IM == 12:
        n_nodes, epsi_nodes, weight_nodes = Monomials_2(N,vcv)

    # 3. Polynomial bases for the test
    X = Ord_Polynomial_N(np.hstack((k,a)),D)

    # 4. Given the solution for capital, compute consumption on the given set of points
    Error =[]
    for i in range(1, P+1):
        # 4.1 Variables in point p
        i = 1
        # N capital stocks of period t
        k0 = k[i-1,:N][np.newaxis,:]
        # N productivity levels of period t
        a0 = a[i-1,:N][np.newaxis,:]
        # Complete (second-degree) polynomial
        X0 = X[i-1, :][np.newaxis,:]

        #4.2 Capital and consumption choices at t

        # Compute a row-vector of capital of period t+1 (chosen at t) using
        # the corresponding capital policy functions
        k1 = X0@bk

        C0 = (A*k0**alpha*a0 - k1 + (1-delta)*k0)@np.ones((N,1))

        c0 = C0@np.ones((1,N))/N
        #4.3 Capital and consumption choices at t+1

        #Compute the next-period productivity levels in each integration node
        a1 = (np.ones((n_nodes,1))@a0)**rho*np.exp(epsi_nodes)

        #Duplicate k1 n_nodes times to create a matrix with n_nodes identical rows
        k1_dupl = np.ones((n_nodes,1))@k1

        # Form a complete polynomial of degree D (at t+1) in the given point 
        X1 = Ord_Polynomial_N(np.hstack((k1_dupl,a1)),D)

        #Compute capital of period t+2 (chosen at t+1) using the second-degree capital policy functions
        k2 = X1@bk

        # Aggregate consumption is computed by summing up individual consumption, which in turn, is found from the individual budget constraints
        C1 = (A*k1_dupl**alpha*a1-k2+ (1-delta)*k1_dupl)@np.ones((N,1))

        c1 = C1@np.ones((1,N))/N

        #5. Approximation errors in point p
        ###################################

        # 5.1 Lagrange multiplier associated with the aggregate resource constraint
        #Compute a country's marginal utility of consumption multiplied by its welfare weight
        MUC0j = np.empty((1,N))
        for j in range(1,N+1):
            MUC0j[0,j-1] = tau*c0[0, j-1]**(-gam)
        # An optimality condition w.r.t. consumption of period t equates 
        # the Lagrange multiplier of the aggregate resource constraint of 
        # period t and each country's marginal utility of consumption 
        # multiplied by its welfare weight; to infer the Lagrange multiplier,  
        # we average across N countries
        lambda0 = np.mean(MUC0j, axis=1)[np.newaxis,:]

        MUC1j = np.empty((n_nodes, N))

        #Compute a country's marginal utility of consumption multiplied by its welfare weight
        for j in range(1,N+1):
            MUC1j[:n_nodes, j-1] = tau*c1[:n_nodes,j-1]**(-gam)
        # Similarly, the Lagrange multiplier of the aggregate resource 
        # constraint of period t+1 is equal to a country's marginal utility 
        # of consumption multiplied by its welfare weight; to infer the 
        # Lagrange multiplier, we average across N countries
        lambda1 = np.mean(MUC1j, axis=1)[:,np.newaxis]

        #5.2 Unit-free Euler-equation errors
        Errors1 = np.zeros((1,N))
        for j in range(1, N+1):
            Errors1[0,j-1] = 1- weight_nodes.conj().T@(beta*lambda1/lambda0*(1-delta+alpha*A*k1[0,j-1]**(alpha-1)*a1[:n_nodes,j-1])[:,np.newaxis])
        # 5.2 Unit-free errors in the optimality conditions w.r.t. consumption
        Errors2 = np.zeros((1,N))
        for j in range(1, N+1):
            Errors2[0,j-1] = 1- lambda0/(tau*c0[0,j-1]**(-gam))
        # 5.3 Unit-free errors in the optimality conditions w.r.t. labor 
        Errors3 = np.zeros((1, N))
        # 5.4 Unit-free approximation error in the aggregate resource constraint
        Errors4 = 1- (c0[0,:N][np.newaxis,:] + k1[0,0:N][np.newaxis,:]-(1-delta)*k0[0,0:N])@np.ones((N,1))/(A*k0[0,0:N][np.newaxis,:]**alpha*a0[0,0:N][np.newaxis,:]@np.ones((N,1)))
        # 5.5 Approximation errors in the capital-accumulation equation
        Errors5 = np.zeros((1,N))

        # combine it into a row
        # First we flatten it in order to store the values in a list
        Error_final = np.hstack((Errors1, Errors2, Errors3, Errors4, Errors5)).flatten()
        Error.append(Error_final)
    #Restore it to 2D
    Error = np.stack(Error, axis=0)

    #6. Mean and maximum approximation errors computed after discarding the "discard" observations
    ##############################################################################################

    # 6.1 Approximation errors across all the optimality conditions
    # Average absolute approximation error 
    Errors_mean = math.log10(np.mean(np.abs(Error[discard:,])))
    # Maximum absolute approximation error
    Errors_max = math.log10(np.max(np.abs(Error[discard:])))
    # 6.2 is not touched as it is not relevent for the results
    end = time.time()
    time_test = end-start

    return Errors_mean, Errors_max, time_test

def GSSA_country(N=2):
    '''
    This function showcase the python implementation from the original matlab codes.
    ------
    Argument:
        N(int):Number of countries. Default is 2

    ------
    Output:
        Time(list): List of run time for each polynominal in Stage 2.
        stage1(float): Time for initial Monte Carol guess (Stage 1).
        Max(list): List of maximum approximation errors for each polynominal.
        Mean(list): List of mean approximation errors for each polynominal.
        Test_time(list): List of time for testing accuracy.
    '''
    #Road map:
    #The 4 main parts follow the construction of the original code in Matlab
    #1: Simulation
    #2: GSSA stage 1:Initial Guess
    #3: GSSA stage 1:Updating Initial Guess
    #4: GSSA stage 2

    #################################
    #Simulation
    #################################

    # User defined parameters
    N = N
    T = 2000
    # Model parameters
    gam = 1
    alpha = 0.36
    beta = 0.99
    delta = 0.025
    rho = 0.95
    sigma = 0.01
    # Variance-covariance matrix
    vcv = sigma**(2)*(np.eye(N)+ np.ones(N))
    # Normalising constant
    A = (1-beta+beta*delta)/alpha/beta
    # welfare weight
    tau = 1
    # The above normalization ensures that steady state
    k = np.ones((T+1,N))
    a = np.ones((1,N))

    #construct productivity levels
    #a20200 = Productivity(T,N,a,sigma,rho)

    # Notice: Similar to the original matlab code, we will be using a prepared dataset by the authors
    # but one can of course uncomment the Productivity and simulate it
    real_path= os.path.join(os.getcwd(), "data\\")
    a20200 = scipy.io.loadmat(real_path+r"aT20200N10.mat").get("a20200")
    a20200 = a20200[:T, :N]

    

    ##########################
    # GSSA stage 1:First Guess
    ##########################
    # Compute a first-degree polynomial solution using the one-node Monte Carlo integration method 
    # (this #solution will be used as an initial guess for the other cases) 

    start = time.time()
    kdamp = 0.1
    dif_1d = 1e+10
    # Initialize the first-degree capital policy functions of N countries 
    bk_1d= np.vstack((np.zeros((1,N)), np.diag(0.9*np.ones(N)),np.diag(0.1*np.ones(N))))
    # Initialize the capital series
    k_old = np.ones((T+1,N))
    # containers
    Y = np.empty((T-1,N))
    # The main iterative cycle of GSSA
    while dif_1d > 1e-4*kdamp:
        # construct starting numpy for value storages
        x = np.empty((T,2*N+1))
        for i in range(T):
            x[i]= np.hstack((1,k[i], a20200[i]))
            k[i+1]= np.hstack((1,k[i], a20200[i]))@bk_1d

        # Compute consumption series 
        C = (A*k[:T]**alpha*a20200[:T]-k[1:T+1]+(1-delta)*k[:T])@np.ones((N,1))
        # Individual consumption is the same for all countries
        c = C@np.ones((1,N))/N
        # Evaluate the percentage (unit-free) difference between the series from the previous and current iterations
        dif_1d = np.mean(abs(1-k/k_old))
        # Monte Carlo realizations of the right side of the Euler equation
        
        for i in range(N):
            Y[:T-1,i]= beta*c[1:T,i]**(-gam)/c[:T-1,i]**(-gam)*(1-delta+alpha*A*k[1:T,i]**(alpha-1)*a20200[1:T,i])*k[1:T,i]
        
        #Compute and update the coefficients of the capital policy functions
        bk_hat_1d = linalg.inv(x[:T-1].conj().T@x[:T-1])@x[:T-1].conj().T@Y[:T-1]
        bk_1d = kdamp*bk_hat_1d+(1-kdamp)*bk_1d
        k_old = k.copy()

    # timing the function
    end = time.time()
    stage1 = end-start

    ######################################
    # GSSA Stage 1: Updating initial guess
    ######################################
    # Compute polynomial solutions of the degrees from one to D_max using one of the following integration methods: 
    # Monte Carlo, Gauss-Hermite product and monomial non-product methods


    # Damping parameter for (fixed-point) iteration on the coefficients of the capital policy functions
    kdamp     = 0.1
    # Set the initial difference between the series from two iterations in the convergence criterion   
    # dif_GSSA_D  = 1e+10

    # The matrix of the polynomial coefficients
    D_max = 5

    # construct container
    npol = np.empty((1,D_max))

    for i in range(1, D_max+1):
        npol[0,i-1]= Ord_Polynomial_N(np.hstack((k[0],a20200[0]))[np.newaxis,], i).shape[1]
    # Choose an integration method here
    ###add more comments here
    IM = 11

    ##################################
    if ((IM>=1) and (IM<=10)):
        n_nodes, epsi_nodes, weight_nodes = GH_Quadrature(IM, N, vcv)
    elif IM == 11:
        n_nodes, epsi_nodes, weight_nodes = Monomials_1(N,vcv)
    elif IM == 12:
        n_nodes, epsi_nodes, weight_nodes = Monomials_2(N,vcv)

    # Choose a regression method and specifications
    #####only 6 methods avaliable here !!!!####

    RM    = 5       
    normalize = 1
    penalty = 7 

    #Compute the polynomial solutions of the degrees from one to D_max
    #Construct empty container for results
    Time = []
    BK = []

    for D in range(1, D_max+1):
        start = time.time()
        # Using the previously computed capital series, compute the initial 
        # guess on the coefficients under the  selected approximation method
        X = Ord_Polynomial_N(np.hstack((k[:T,:],a20200[:T,:])),D)
        bk_D = Num_Stab_Approx(X[:T-1,:],Y[:T-1,:], RM, penalty, normalize)
        k_old = np.ones((T+1,N))
        dif_GSSA_D  = 1e+10
        while dif_GSSA_D > 1e-4/10**D*kdamp:
            #for checking purposes
            #print(dif_GSSA_D)

            #Generate time series of capital
            for i in range(1, T+1):
                X[i-1,] = Ord_Polynomial_N(np.hstack((k[i-1,:], a20200[i-1, :]))[np.newaxis,], D)
                k[i,] = X[i-1,][np.newaxis,]@bk_D
            
            X = np.nan_to_num(X)
            #15.2.2Compute consumption series of all countries:
            ###################################################
            # N current capital stocks  
            k0 = k[0:T,]
            # N current productivity levels 
            a0 = a20200[0:T,]
            # N next-period capital stocks
            k1 = k[1:T+1,]
            # Aggregate consumption is computed by summing up individual consumption, which in turn, is found from the individual budget constraints
            C = (A*k0**alpha*a0 - k1+ (1-delta)*k0)@np.ones((N,1))

            # Individual consumption is the same for all countries, check JMM (2011)
            c = C@np.ones((1,N))/N

            # Approximate the conditional expectations for t=1,...T-1 using the integration method chosen

            # The one-node Monte Carlo integration method approximates the 
            # values of the conditional expectations, Y, in the Euler equation with
            # the realization of the integrand in the next period
            if IM == 0:
                for i in range(1, N+1):
                    Y[0:T-1,i-1][:,np.newaxis] = beta*c[1:T,i-1][:,np.newaxis]**(-gam)/c[0:T-1, i-1][:,np.newaxis]**(-gam)*(1-delta+alpha*A*k[1:T,i-1][:,np.newaxis]**(alpha-1)*a20200[1:T,i-1][:,np.newaxis])*k[1:T,i-1][:,np.newaxis]

                Y = np.nan_to_num(Y)
            # Deterministic integration methods approximate the values of 
            # conditional expectations, Y, in the Euler equation as a weighted average 
            # of the values of the integrand in the given nodes with the given weights 
            else:
                Y = np.zeros((T,N))
                for i in range(1, n_nodes+1):
                    # Compute the next-period productivity levels for each integration node using condition (C3) in JMM (2011)
                    a1 = a20200[0:T,:]**rho*np.exp(np.ones((T,1))@epsi_nodes[i-1,:][np.newaxis,])

                    # Compute capital of period t+2 (chosen at t+1) using the capital policy functions
                    k2 = Ord_Polynomial_N(np.hstack((k1, a1)),D)@bk_D

                    # C is computed by summing up individual consumption, which in turn, is from the individual budget constraints
                    C1 = (A*k1**alpha*a1-k2+(1-delta)*k1)@np.ones((N,1))

                    # Compute next-period individual consumption for N countries 
                    c1 = C1@np.ones((1,N))/N
                    for i in range(1, N+1):
                        Y[0:T,i-1][:,np.newaxis]= Y[0:T,i-1][:,np.newaxis]+ weight_nodes[i-1, 0]*beta*c1[0:T, i-1][:,np.newaxis]**(-gam)/c[0:T,i-1][:,np.newaxis]**(-gam)*(1-delta+alpha*A*k1[0:T,i-1][:,np.newaxis]**(alpha-1)*a1[0:T,i-1][:,np.newaxis])*k1[0:T, i-1][:,np.newaxis]

                Y = np.nan_to_num(Y)
            # Evaluate the percentage (unit-free) difference between the 
            # capital series from the previous and current iterations
            dif_GSSA_D = np.mean(abs(1-k/k_old))

            # 15.2.5 Compute and update the coefficients of the capital policy 
            # functions 
 
            # Compute new coefficients of the capital 
            # policy functions using the chosen 
            # approximation method
            bk_hat_D = Num_Stab_Approx(X[0:T-1,:], Y[0:T-1,:], RM, penalty, normalize)
            bk_D = kdamp*bk_hat_D+ (1-kdamp)*bk_D

            #15.2.6 Store the capital series 
            #--------------------------------
            k_old = k.copy()

        # The GSSA output for the polynomial solution of degree D
        end = time.time()
        BK.append(bk_D)
        Time.append(end-start)

    ##############
    # GSSA stage 2
    ##############
    # Choose the simulation length for the test on a stochastic simulation, T_test<=10,200
    T_test = 10200
    a20200 = scipy.io.loadmat(real_path+r"aT20200N10.mat").get("a20200")

    # Restrict the series of the productivity levels for testing
    a_test = a20200[T:T+T_test,0:N]

    # Choose an integration method for evaluating accuracy of solutions
    IM_test = 11

    # Compute errors on a stochastic simulation for the GSSA polynomial solution of degrees D=1,...,D_max
    Max = []
    Mean = []
    Test_time = []
    for D in range(1, D_max+1):
        # Simulate the time series solution under the given capital-
        # policy-function coefficients, BK(:,:,D) with D=1,...,D_max 
        bk = BK[D-1]
        # Initial condition for capital (equal to steady state)
        k_test = np.ones((1,N))

        for j in range(1, T_test):
            X_test = Ord_Polynomial_N(np.hstack((k_test[j-1][np.newaxis,],a_test[j-1,][np.newaxis,])), D)
            k_test = np.vstack((k_test, X_test@bk))
    
        # Errors across 10,000 points on a stochastic simulation
        discard = 200
        Errors_mean, Errors_max, time_test = Accuracy_Test_N(k_test, a_test, bk, D, IM_test, alpha, gam, delta, beta, A, tau, rho, vcv, discard)
        Max.append(Errors_max)
        Mean.append(Errors_mean)
        Test_time.append(time_test)

    return Time, stage1, Max, Mean, Test_time

def GSSA_country_hd(N=2, D_max=1):
    '''
    This function for high dimension calculation.
    ------
    Argument:
        N(int):Number of countries. Default is 2
        D_max(int):Number of maximum polynomial degree. Default is 1

    ------
    Output:
        Time(list): List of run time for each polynominal in Stage 2.
        stage1(float): Time for initial Monte Carol guess (Stage 1).
        Max(list): List of maximum approximation errors for each polynominal.
        Mean(list): List of mean approximation errors for each polynominal.
        Test_time(list): List of time for testing accuracy.
    '''
    # All comments are similar following the original code


    # User defined parameters
    N = N
    T = 2000
    # Model parameters
    gam = 1
    alpha = 0.36
    beta = 0.99
    delta = 0.025
    rho = 0.95
    sigma = 0.01
    # Variance-covariance matrix
    vcv = sigma**(2)*(np.eye(N)+ np.ones(N))
    # Normalising constant
    A = (1-beta+beta*delta)/alpha/beta
    # welfare weight
    tau = 1
    #The above normalization ensures that steady state
    k = np.ones((T+1,N))
    a = np.ones((1,N))


    if N<=10:
        real_path= os.path.join(os.getcwd(), "data\\")
        a20200 = scipy.io.loadmat(real_path+r"aT20200N10.mat").get("a20200")
        a20200 = a20200[:T, :N]
    else:
        #construct productivity levels
        a20200 = Productivity(T,N,a,sigma,rho)
    
    ############################
    #GSSA stage 1: Initial guess
    ############################

    start = time.time()
    kdamp = 0.1
    dif_1d = 1e+10
    # Initialize the first-degree capital policy functions of N countries 
    bk_1d= np.vstack((np.zeros((1,N)), np.diag(0.9*np.ones(N)),np.diag(0.1*np.ones(N))))
    # Initialize the capital series
    k_old = np.ones((T+1,N))
    # containers
    Y = np.empty((T-1,N))
    # The main iterative cycle of GSSA
    while dif_1d > 1e-4*kdamp:
        #construct starting numpy for value storages
        x = np.empty((T,2*N+1))
        for i in range(T):
            x[i]= np.hstack((1,k[i], a20200[i]))
            k[i+1]= np.hstack((1,k[i], a20200[i]))@bk_1d

        #Compute consumption series 
        C = (A*k[:T]**alpha*a20200[:T]-k[1:T+1]+(1-delta)*k[:T])@np.ones((N,1))
        #Individual consumption is the same for all countries
        c = C@np.ones((1,N))/N
        #Evaluate the percentage (unit-free) difference between the series from the previous and current iterations
        dif_1d = np.mean(abs(1-k/k_old))
        #Monte Carlo realizations of the right side of the Euler equation
        
        for i in range(N):
            Y[:T-1,i]= beta*c[1:T,i]**(-gam)/c[:T-1,i]**(-gam)*(1-delta+alpha*A*k[1:T,i]**(alpha-1)*a20200[1:T,i])*k[1:T,i]
        
        #Compute and update the coefficients of the capital policy functions
        bk_hat_1d = linalg.inv(x[:T-1].conj().T@x[:T-1])@x[:T-1].conj().T@Y[:T-1]
        bk_1d = kdamp*bk_hat_1d+(1-kdamp)*bk_1d
        k_old = k.copy()

    #timing the function
    end = time.time()
    stage1 = end-start

    #################################
    #GSSA Stage 1: Updating the guess
    ################################

    # Damping parameter for (fixed-point) iteration on the coefficients of the capital policy functions
    kdamp     = 0.1
    #Set the initial difference between the series from two iterations in the convergence criterion   
    #dif_GSSA_D  = 1e+10

    # The matrix of the polynomial coefficients
    D_max = D_max

    #construct container
    npol = np.empty((1,D_max))

    for i in range(1, D_max+1):
        npol[0,i-1]= Ord_Polynomial_N(np.hstack((k[0],a20200[0]))[np.newaxis,], i).shape[1]
    # Choose an integration method here
    ###add more comments here
    IM = 11

    ##################################
    if ((IM>=1) and (IM<=10)):
        n_nodes, epsi_nodes, weight_nodes = GH_Quadrature(IM, N, vcv)
    elif IM == 11:
        n_nodes, epsi_nodes, weight_nodes = Monomials_1(N,vcv)
    elif IM == 12:
        n_nodes, epsi_nodes, weight_nodes = Monomials_2(N,vcv)

    #14. Choose a regression method and specifications
    #####only 6 methods avaliable here !!!!####
    RM    = 5       
    normalize = 1
    penalty = 7 

    
    #Construct empty container for results
    Time = []
    BK = []
    #Compute the polynomial solutions of the degrees from one to D_max
    for D in range(1, D_max+1):
        start = time.time()
        # Using the previously computed capital series, compute the initial 
        # guess on the coefficients under the  selected approximation method
        X = Ord_Polynomial_N(np.hstack((k[:T,:],a20200[:T,:])),D)
        bk_D = Num_Stab_Approx(X[:T-1,:],Y[:T-1,:], RM, penalty, normalize)
        k_old = np.ones((T+1,N))
        dif_GSSA_D  = 1e+10
        while dif_GSSA_D > 1e-4/10**D*kdamp:
            #for checking purposes
            #print(dif_GSSA_D)

            #Generate time series of capital
            for i in range(1, T+1):
                X[i-1,] = Ord_Polynomial_N(np.hstack((k[i-1,:], a20200[i-1, :]))[np.newaxis,], D)
                k[i,] = X[i-1,][np.newaxis,]@bk_D
            
            X = np.nan_to_num(X)
            # Compute consumption series of all countries:
            ###################################################
            # N current capital stocks  
            k0 = k[0:T,]
            # N current productivity levels 
            a0 = a20200[0:T,]
            # N next-period capital stocks
            k1 = k[1:T+1,]
            # Aggregate consumption is computed by summing up individual consumption, which in turn, is found from the individual budget constraints
            C = (A*k0**alpha*a0 - k1+ (1-delta)*k0)@np.ones((N,1))

            # Individual consumption is the same for all countries, check JMM (2011)
            c = C@np.ones((1,N))/N

            # Approximate the conditional expectations for t=1,...T-1 using the integration method chosen

            # The one-node Monte Carlo integration method approximates the 
            # values of the conditional expectations, Y, in the Euler equation with
            # the realization of the integrand in the next period
            if IM == 0:
                for i in range(1, N+1):
                    Y[0:T-1,i-1][:,np.newaxis] = beta*c[1:T,i-1][:,np.newaxis]**(-gam)/c[0:T-1, i-1][:,np.newaxis]**(-gam)*(1-delta+alpha*A*k[1:T,i-1][:,np.newaxis]**(alpha-1)*a20200[1:T,i-1][:,np.newaxis])*k[1:T,i-1][:,np.newaxis]

                Y = np.nan_to_num(Y)
            #  Deterministic integration methods approximate the values of 
            # conditional expectations, Y, in the Euler equation as a weighted average 
            # of the values of the integrand in the given nodes with the given weights 
            else:
                Y = np.zeros((T,N))
                for i in range(1, n_nodes+1):
                    # Compute the next-period productivity levels for each integration node using condition (C3) in JMM (2011)
                    a1 = a20200[0:T,:]**rho*np.exp(np.ones((T,1))@epsi_nodes[i-1,:][np.newaxis,])

                    # Compute capital of period t+2 (chosen at t+1) using the capital policy functions
                    k2 = Ord_Polynomial_N(np.hstack((k1, a1)),D)@bk_D

                    # C is computed by summing up individual consumption, which in turn, is from the individual budget constraints
                    C1 = (A*k1**alpha*a1-k2+(1-delta)*k1)@np.ones((N,1))

                    # Compute next-period individual consumption for N countries 
                    c1 = C1@np.ones((1,N))/N
                    for i in range(1, N+1):
                        Y[0:T,i-1][:,np.newaxis]= Y[0:T,i-1][:,np.newaxis]+ weight_nodes[i-1, 0]*beta*c1[0:T, i-1][:,np.newaxis]**(-gam)/c[0:T,i-1][:,np.newaxis]**(-gam)*(1-delta+alpha*A*k1[0:T,i-1][:,np.newaxis]**(alpha-1)*a1[0:T,i-1][:,np.newaxis])*k1[0:T, i-1][:,np.newaxis]

                Y = np.nan_to_num(Y)
            # Evaluate the percentage (unit-free) difference between the 
            # capital series from the previous and current iterations
            #-------------------------------------------------------
            dif_GSSA_D = np.mean(abs(1-k/k_old))

            # Compute and update the coefficients of the capital policy functions 

            # Compute new coefficients of the capital 
            # policy functions using the chosen 
            # approximation method
            bk_hat_D = Num_Stab_Approx(X[0:T-1,:], Y[0:T-1,:], RM, penalty, normalize)
            bk_D = kdamp*bk_hat_D+ (1-kdamp)*bk_D

            # Store the capital series 
            k_old = k.copy()

        # The GSSA output for the polynomial solution of degree D
        end = time.time()
        BK.append(bk_D)
        Time.append(end-start)

    #############
    #GSSA stage 2
    #############
    # Choose the simulation length for the test on a stochastic simulation, T_test<=10,200
    T_test = 10200
    if N<=10:
        real_path= os.path.join(os.getcwd(), "data\\")
        a20200 = scipy.io.loadmat(real_path+r"aT20200N10.mat").get("a20200")
        a20200 = a20200[:T, :N]
    else:
        #construct productivity levels
        a20200 = Productivity(T_test,N,a,sigma,rho)

    # Restrict the series of the productivity levels for the test
    a_test = a20200[:T_test,:N]

    # Choose an integration method for evaluating accuracy of solutions
    IM_test = 11

    # Compute errors on a stochastic simulation for the GSSA polynomial solution of degrees D=1,...,D_max
    Max = []
    Mean = []
    Test_time = []
    for D in range(1, D_max+1):
        # Simulate the time series solution under the given capital-
        # policy-function coefficients, BK(:,:,D) with D=1,...,D_max 
        bk = BK[D-1]
        # Initial condition for capital (equal to steady state)
        k_test = np.ones((1,N))

        for j in range(1, T_test):
            X_test = Ord_Polynomial_N(np.hstack((k_test[j-1][np.newaxis,],a_test[j-1,][np.newaxis,])), D)
            k_test = np.vstack((k_test, X_test@bk))
    
        # Errors across 10,000 points on a stochastic simulation
        discard = 200
        Errors_mean, Errors_max, time_test = Accuracy_Test_N(k_test, a_test, bk, D, IM_test, alpha, gam, delta, beta, A, tau, rho, vcv, discard)
        Max.append(Errors_max)
        Mean.append(Errors_mean)
        Test_time.append(time_test)

    return Time, stage1, Max, Mean, Test_time


def GSSA_country_df(N=2, Cache=True):
    '''
    This function prepares a dataframe for result showcasing.
    ------
    Argument:
        N(int): Number of countries. Default is 2
        Cache(boolean): Cache the result or not. Default is True.
    ------
    Output:
        df(pandas DF): the result of 
    '''
    if ((N==1)&(Cache==True)):
        real_path= os.path.join(os.getcwd(), "cache\\")
        df = pd.read_csv(real_path+ r"countries_1.csv")
    elif ((N==2)&(Cache==True)):
        real_path= os.path.join(os.getcwd(), "cache\\")
        df = pd.read_csv(real_path+ r"countries_2.csv")
    
    elif N<=2:
        Time, stage1, Max, Mean, Test_time = GSSA_country(N=N)
        # DF for the results
        df = pd.DataFrame({"Polynomial Degree":[i for i in range(1,6)],
            "Mean Error":Mean,
            "Max Error":Max,
            "Total Time":[Time[i]+stage1+Test_time[i] for i in range(5)],
            "Number of countries":str(N)
            })
    elif ((N>2)&(N<=10)):
        Time, stage1, Max, Mean, Test_time = GSSA_country_hd(N=N,D_max=2)
        # DF for the results
        df = pd.DataFrame({"Polynomial Degree":[i for i in range(1,3)],
            "Mean Error":Mean,
            "Max Error":Max,
            "Total Time":[Time[i]+stage1+Test_time[i] for i in range(2)],
            "Number of countries":str(N)
            })

    else:
        Time, stage1, Max, Mean, Test_time = GSSA_country_hd(N=N,D_max=1)
        # DF for the results
        df = pd.DataFrame({"Polynomial Degree":[i for i in range(1,2)],
            "Mean Error":Mean,
            "Max Error":Max,
            "Total Time":[Time[i]+stage1+Test_time[i] for i in range(1)],
            "Number of countries":str(N)
            })

    return df
