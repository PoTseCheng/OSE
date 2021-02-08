import numpy as np
import pandas as pd
import math
import time
import scipy.linalg
from numpy.linalg import inv
from numpy.linalg import svd
from scipy.optimize import linprog
from scipy.special import erfinv
from numpy.random import rand
from numpy import sqrt

##########################################
#General information regarding this file
##########################################
#Remark 1:
#This is the python implementation for the paper "Numerically Stable and Accurate Stochastic Simulation Approaches
#for Solving Dynamic Economic Models" by Kenneth L. Judd, Lilia Maliar and Serguei Maliar, (2011), henceforth JMM (2011)

#Remark 2:
#This file contains codes from my initial attemps of translating matlab codes from the original authors. 
#These codes are therefore less stylish compared to the codes in country.py which were done in later stages.


def randn2(*args,**kwargs):
    '''
    The randn function that matlab uses.
    ------
    '''
    uniform = rand(*args, **kwargs)
    return sqrt(2) * erfinv(2 * uniform - 1)
 

def reading(x):
    '''
    Homemade reader for csv file produced by homemade converter.
    --------
    Parameters:
        x(str):Csv file from converter

    Returns:
        df(Pandas dataframe)
    '''
    df=pd.read_csv(x)
    add=df.columns.tolist()
    df.loc[-1] = add 
    df.index = df.index + 1  
    df = df.sort_index()
    df.columns=[x.replace(".csv","")]

    return df

def Ord_Herm_Pol_1(z,D,PF,zb):
    '''
    The function constructs the basis functions of complete ordinary and Hermite polynomials of the degrees from one to five for the two-dimensional (two-state variables) cases.
    ----------
    Arguments:
        z(numpy array): Data points on which the polynomial basis functions must be constructed.
        D(int): The degree of the polynomial whose basis functions must be constructed. Can be 1~5
        PF(binary): The polynomial family chosen: 0=Ordinary, 1=Hermite.
        zb(numpy array): Matrix of means and standard deviations of state variables. 
    
    --------
    Output:
        basis(numpy array): The array of basis functions of a complete polynomial of the given degree.
    '''
    
    num_rows = z.shape[0]
    if PF == 1:
        zc1 = (z[:, 0]-zb[0, 0])/zb[1, 0]
        zc2 = (z[:, 1]-zb[0, 1])/zb[1, 1]
        # p1,...,p5 are the vectors obtained by evaluating the Hermite
        p1 = zc1
        p2 = np.power(zc1, 2)-1
        p3 = np.power(zc1, 3)-3*zc1
        p4 = np.power(zc1, 4)-6*(np.power(zc1, 2))+3
        p5 = np.power(zc1, 5)-10*(np.power(zc1, 3))+15*zc1
        # q1,...,q5 are the vectors obtained by evaluating the Hermite 
        q1 = zc2
        q2 = np.power(zc2, 2)-1
        q3 = np.power(zc2, 3)-3*zc2
        q4 = np.power(zc2, 4)-6*(np.power(zc2, 2))+3
        q5 = np.power(zc2, 5)-10*(np.power(zc2, 3))+15*zc2
        if D == 1:
            basis = np.concatenate((np.ones((num_rows, 1)),
            p1.reshape(num_rows, 1),
            q1.reshape(num_rows, 1)
            ), axis=1)
            return basis

                
        elif D == 2:
            basis = np.concatenate((np.ones((num_rows, 1)),
            p1.reshape(num_rows, 1),
            q1.reshape(num_rows, 1),
            p2.reshape(num_rows, 1),
            np.multiply(p1, q1).reshape(num_rows, 1),
            q2.reshape(num_rows, 1)
            ), axis=1)
            return basis


        elif D == 3:
            basis = np.concatenate((np.ones((num_rows, 1)),
            p1.reshape(num_rows, 1),
            q1.reshape(num_rows, 1),
            p2.reshape(num_rows, 1),
            np.multiply(p1, q1).reshape(num_rows, 1),
            q2.reshape(num_rows, 1),
            p3.reshape(num_rows, 1),
            np.multiply(p2, q1).reshape(num_rows, 1),
            np.multiply(p1, q2).reshape(num_rows, 1),
            q3.reshape(num_rows, 1)
            ), axis=1)
            return basis


        elif D == 4:
            basis = np.concatenate((np.ones((num_rows, 1)),
            p1.reshape(num_rows, 1),
            q1.reshape(num_rows, 1),
            p2.reshape(num_rows, 1),
            np.multiply(p1, q1).reshape(num_rows, 1),
            q2.reshape(num_rows, 1),
            p3.reshape(num_rows, 1),
            np.multiply(p2, q1).reshape(num_rows, 1),
            np.multiply(p1, q2).reshape(num_rows, 1),
            q3.reshape(num_rows, 1),
            p4.reshape(num_rows, 1),
            np.multiply(p3, q1).reshape(num_rows, 1),
            np.multiply(p2, q2).reshape(num_rows, 1),
            np.multiply(p1, q3).reshape(num_rows, 1),
            q4.reshape(num_rows, 1)
            ), axis=1)
            return basis


        elif D == 5:
            basis = np.concatenate((np.ones((num_rows, 1)),
            p1.reshape(num_rows, 1),
            q1.reshape(num_rows, 1),
            p2.reshape(num_rows, 1),
            np.multiply(p1, q1).reshape(num_rows, 1),
            q2.reshape(num_rows, 1),
            p3.reshape(num_rows, 1),
            np.multiply(p2, q1).reshape(num_rows, 1),
            np.multiply(p1, q2).reshape(num_rows, 1),
            q3.reshape(num_rows, 1),
            p4.reshape(num_rows, 1),
            np.multiply(p3, q1).reshape(num_rows, 1),
            np.multiply(p2, q2).reshape(num_rows, 1),
            np.multiply(p1, q3).reshape(num_rows, 1),
            q4.reshape(num_rows, 1),
            p5.reshape(num_rows, 1),
            np.multiply(p4, q1).reshape(num_rows, 1),
            np.multiply(p3, q2).reshape(num_rows, 1),
            np.multiply(p2, q3).reshape(num_rows, 1),
            np.multiply(p1, q4).reshape(num_rows, 1),
            q5.reshape(num_rows, 1)
            ), axis=1)
            return basis
            
        # If the polynomial family chosen is ordinary##
    else:
        zc1 = z[:, 0] # No normalization
        zc2 = z[:, 1] # No normalization
        #p1,...,p5 are the vectors obtained by evaluating the ordinary
        p1 = zc1
        p2 = np.power(zc1, 2)
        p3 = np.power(zc1, 3)
        p4 = np.power(zc1, 4)
        p5 = np.power(zc1, 5)
        #q1,...,q5 are the vectors obtained by evaluating the ordinary
        q1 = zc2
        q2 = np.power(zc2, 2)
        q3 = np.power(zc2, 3)
        q4 = np.power(zc2, 4)
        q5 = np.power(zc2, 5)
        # Construct the matrix of the basis functions
        if D == 1:
            basis = np.concatenate((np.ones((num_rows, 1)),
            p1.reshape(num_rows, 1),
            q1.reshape(num_rows, 1)
            ), axis=1)

            return basis

        elif D == 2:
            basis = np.concatenate((np.ones((num_rows, 1)),
            p1.reshape(num_rows, 1),
            q1.reshape(num_rows, 1),
            p2.reshape(num_rows, 1),
            np.multiply(p1, q1).reshape(num_rows, 1),
            q2.reshape(num_rows, 1)
            ), axis=1)

            return basis

        elif D == 3:
            basis = np.concatenate((np.ones((num_rows, 1)),
            p1.reshape(num_rows, 1),
            q1.reshape(num_rows, 1),
            p2.reshape(num_rows, 1),
            np.multiply(p1, q1).reshape(num_rows, 1),
            q2.reshape(num_rows, 1),
            p3.reshape(num_rows, 1),
            np.multiply(p2, q1).reshape(num_rows, 1),
            np.multiply(p1, q2).reshape(num_rows, 1),
            q3.reshape(num_rows, 1)
            ), axis=1)

            return basis
        elif D == 4:
            basis = np.concatenate((np.ones((num_rows, 1)),
            p1.reshape(num_rows, 1),
            q1.reshape(num_rows, 1),
            p2.reshape(num_rows, 1),
            np.multiply(p1, q1).reshape(num_rows, 1),
            q2.reshape(num_rows, 1),
            p3.reshape(num_rows, 1),
            np.multiply(p2, q1).reshape(num_rows, 1),
            np.multiply(p1, q2).reshape(num_rows, 1),
            q3.reshape(num_rows, 1),
            p4.reshape(num_rows, 1),
            np.multiply(p3, q1).reshape(num_rows, 1),
            np.multiply(p2, q2).reshape(num_rows, 1),
            np.multiply(p1, q3).reshape(num_rows, 1),
            q4.reshape(num_rows, 1)
            ), axis=1)

            return basis
        elif D == 5:
            
            basis = np.concatenate((np.ones((num_rows, 1)),
            p1.reshape(num_rows, 1),
            q1.reshape(num_rows, 1),
            p2.reshape(num_rows, 1),
            np.multiply(p1, q1).reshape(num_rows, 1),
            q2.reshape(num_rows, 1),
            p3.reshape(num_rows, 1),
            np.multiply(p2, q1).reshape(num_rows, 1),
            np.multiply(p1, q2).reshape(num_rows, 1),
            q3.reshape(num_rows, 1),
            p4.reshape(num_rows, 1),
            np.multiply(p3, q1).reshape(num_rows, 1),
            np.multiply(p2, q2).reshape(num_rows, 1),
            np.multiply(p1, q3).reshape(num_rows, 1),
            q4.reshape(num_rows, 1),
            p5.reshape(num_rows, 1),
            np.multiply(p4, q1).reshape(num_rows, 1),
            np.multiply(p3, q2).reshape(num_rows, 1),
            np.multiply(p2, q3).reshape(num_rows, 1),
            np.multiply(p1, q4).reshape(num_rows, 1),
            q5.reshape(num_rows, 1)
            ), axis=1)
            return basis

def Num_Stab_Approx(x,y,RM,penalty,normalised):
    '''
    This function implements the approximation methods mentioned in Judd et al. (2011)
    
    --------
    Arguments:
        x(2D numpy array): Matrix of dependent variables in a regression.
        y(2D numpy array): Matrix of independent variables.
        RM(int): Regression (approximation) method, from 1 to 7. RLAD-DP not included.
        penalty(int): Regularisation parameter for a regularisation methods.
        normalised(int): Optional parameter to normalise the data or not. 1 or 0.

    ---------
    Outputs:
        B(2D numpy array): Matrix of the regression coefficients.

    ---------
    Additional info:
        Here, LAD-DP is included, however this only served as a showcase in the notebook.
    
    '''

    # Step 1: Compute the dimensionality of the data
    #-------------------------------------------------
    T = x.shape[0]
    n = x.shape[1]
    N = y.shape[1] 

    # Step 2: Normalisation
    #----------------------------------
    if ((normalised == 1) or (RM >= 5)):
        X1 = np.divide(
            (x[:,1:n] - np.matmul(np.ones((T, 1)), x[:,1:n].mean(axis=0).reshape(1, n-1))),
            np.matmul(np.ones((T, 1)), np.array([np.std(x[:,i], ddof=1) for i in range(1,n)]).reshape(1,n-1))
        )
        X1 = X1.astype(float)
    
        Y1 = np.divide(
            (y - np.ones((T,1))*np.mean(y)),
            np.matmul(np.ones((T,1)), np.std(y, ddof=1).reshape(1,1))
            )
        Y1 = Y1.astype(float)
        n1 = n-1 # Number of coefficients in a regression with normalized data is reduced by 1 (no intercept)
    else:
        #Leave the values without any changes
        X1 = x
        Y1 = y         
        n1 = n
    
    # Step 3: choosing an approximate method from the paper
    #------------------------------------------------------

    #simple OLS
    if RM == 1:
        B = np.matmul(
            inv(np.matmul(X1.transpose(),X1).astype(float)), #need to convert the output to float or else return error
            np.matmul(X1.transpose(),Y1)
            )

    #LS-SVD
    elif RM == 2:
        U, S, Vh = svd(X1, full_matrices=False)
        V = Vh.T
        S_inv = np.diag(1/S)
        B = np.matmul(np.matmul(np.matmul(V,S_inv),U.transpose()),Y1)
    
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
            #do not change the method here, or else will result in memory overload
            B.append(list(result.x[0:n1]))
        B = np.asarray(B).T
    
    #LAD-DP (This method is unavaliable)
    elif RM == 4:
        BND = [(-1, 1)]*T
        #, bounds= BND
        f = -Y1[:,0]
        #specify the Aeq an beq
        Aeq = X1.conj().transpose()
        B =[]
        #solve the equation
        beq = np.asarray([0]*n1)
        result = linprog(f, A_eq = Aeq, b_eq = beq, bounds= BND, method="highs-ipm")
        B.append(result.x[0:n1])
    #####Despite we can get an answer from Scipy, this is not the answer we want, as we want the lagrange multiplier for coefficient
    
    #RLS-Tikhonov
    elif RM == 5:
        B = np.matmul(
            inv(
                np.matmul(X1.transpose(),X1) + (T/n1)*np.identity(n1)*10**(penalty)
               ),
            np.matmul(X1.transpose(),Y1)
            )

    # RLS-TSVD
    elif RM == 6:
        U, S, Vh = svd(X1, full_matrices=False)
        V = Vh.T
        r = np.count_nonzero(np.divide(np.diag(S).max(), np.diag(S))<= 10**(penalty))
        Sr_inv = np.zeros((n1,n1))
        Sr_inv[0:r, 0:r]= np.diag(np.divide(1., S[:r]))
        B = V@Sr_inv@U.conj().T@Y1

    #LADPP
    elif RM == 7:
        #we can just use the default setting from scipy as the lower and upper will be the same
        f= np.vstack((10**penalty*np.ones((n1*2,1))*T/n1, np.ones((2*T,1))))
        Aeq= np.c_[X1,-X1, np.eye(T), -np.eye(T)]
        B =[]
        #solve the equation
        for i in range(N):
            beq = Y1[:,i]
            result = linprog(f, A_eq = Aeq, b_eq = beq, method="highs-ipm")
            B.append(list(result.x[0:n1]-result.x[n1:2*n1]))
        B = np.asarray(B).T

    #RM == 8 is unavaliable(see notebook)

    #Step 4: Infer the regression coefficients in the original regression with unnormalised data
    #-----------------------------------------------------------------------------------------

    if ((normalised == 1) or (RM >= 5)):
        B2 = np.multiply((1/np.asarray([np.std(x[:,i], ddof=1) for i in range(1,n)])).reshape((n1,1))*np.std(y, ddof=1),B)
        B1 = (y.mean()- np.matmul(np.asarray([x[:,i].mean() for i in range(1,n)]),B2)).reshape(1,1)
        B = np.concatenate((B1,B2))
    #The codes are improved greatly in the implementation, see country.py
    
    return B


def GH_Quadrature(Qn, N, vcv):
    '''
    This function constructs integration nodes and weights under Gauss-Hermite quadrature (product) integration rule with Qn<=10 nodes in each of N dimensions
    ---------
    Arguments:
        Qn(int): The number of nodes in each dimension, notice that Qn must be between 1 to 10.
        N(int): The number of random variables.
        vcv(2D numpy array): The predefined covariance matrix.

    Outputs:
        n_nodes(int):Total number of integration nodes.
        epsi_nodes(2D numpy array):Integration nodes.
        weight_nodes(2D numpy array):Integration weights.
    
    '''

    # For the following codes:
    # Qn : Number of nodes in each dimension; Qn <=10
    # eps: Set of integration nodes
    # weight: Set of integration weights    
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
    else: # The default option is given Qn=10
        Qn = 10 
        eps = np.array([3.436159118837738,2.532731674232790,1.756683649299882,1.036610829789514,0.3429013272237046,-0.3429013272237046,-1.036610829789514,-1.756683649299882,-2.532731674232790,-3.436159118837738])
        weight = np.array([7.640432855232621e-06,0.001343645746781233,0.03387439445548106,0.2401386110823147,0.6108626337353258,0.6108626337353258,0.2401386110823147,0.03387439445548106,0.001343645746781233,7.640432855232621e-06])
    
    # N-dimensional integration nodes and weights for N uncorrelated normally distributed random variables with zero mean and unit variance
    n_nodes = Qn**N

    # Construct containers
    z1 = np.ones([n_nodes,N]).astype(float)
    w1i = np.ones([n_nodes,N]).astype(float)
    w1 = np.ones([n_nodes,1]).astype(float)
    for i in range(N):
        z1[:,i]=np.tile(np.repeat(eps,10**(i)), 10**(N-i-1))
        w1i[:,i]=np.tile(np.repeat(weight,10**(i)), 10**(N-i-1))
    
    for i in range(N):
        w1[:,0]*=w1i[:,i]
    
    # Integration nodes preparations
    z= math.sqrt(2)*z1
    # Integration weights
    weight_nodes= w1/(math.sqrt(math.pi)**N)
    # Integration nodes for condition (B.6) in paper 
    epsi_nodes = np.matmul(z,scipy.linalg.cholesky(vcv))

    return n_nodes,epsi_nodes,weight_nodes


def Accuracy_Test_1(sigma,rho,beta,gam,alpha,delta,k,a,bk,D,IM,PF,zb,discard):
    '''
    '''
    start = time.time()
    n_nodes,epsi_nodes, weight_nodes = GH_Quadrature(Qn=IM, N=1, vcv=sigma**2)
    # Initialise empty list
    errors = []
    for i in range(len(a)):
        k0 = float(k[i])
        a0 = float(a[i])
        k1 = float(Ord_Herm_Pol_1(np.array([k0,a0]).reshape((1,2)), D, PF=PF, zb=zb)@bk)
        c0 = k0**alpha*a0 - k1+ (1-delta)*k0
        a1 = a0**rho*np.exp(epsi_nodes)
        k1_dupl1 = k1*np.ones((n_nodes,1))
        x1 = Ord_Herm_Pol_1(np.concatenate((k1_dupl1, a1), axis=1),D,PF,zb)
        k2 = np.matmul(x1, bk)
        c1 = k1**alpha*a1 - k2 + (1-delta)*k1
        errors.append(weight_nodes.conj().transpose()@(beta*c1**(-gam)/c0**(-gam)*(1-delta+alpha*a1*k1**(alpha-1))) -1) #testing new approach
    Error_means = math.log10(np.mean(np.mean(np.abs(errors[discard:]))))
    Error_max = math.log10(np.max(np.max(np.abs(errors[discard:]))))
    end = time.time()
    time_test = end - start
    return Error_means, Error_max, time_test

def GSSA_main_cycle(T, gam, alpha, beta, delta, kdamp, dif_GSSA_1d, a, bk_1d, k_old, k, checker = 0):
    '''
    Stage 1 of 1 agent model
    --------------
    Parameters:
        T (int): Simulation point
        gam (float): Utility-function parameter.
        alpha(float): Capital share in output.
        beta(float): Discount factor.
        delta(float): Depreciation rate.
        kdamp(float): Damping parameter for (fixed-point) iteration on the coefficients of the capital policy function.
        dif_GSSA_1d(int): Initial difference between the series from two iterations in the convergence criterion.
        a(numpy array): Initial condition for the productivity level.
        bk_1d(numpy array): The coefficients of the capital policy function.
        k_old(numpy array): The series of next-period capital.
        k(numpy array): Initial condition for capital.
        checker (binary): A checker for the while loop, used to validate if the loop is running correctly or not.

    Returns:
        y(numpy array): Monte Carlo realizations of the Euler equation.
    '''
    # The main iterative cycle of GSSA

    # Checking convergence
    while dif_GSSA_1d > 1e-4*kdamp:
        x=[]
        # Generate time series of capital
        for i in range(0,T):
            # The basis functions of the first-degree polynomial at time t
            x.append([1,k[i],a[i]])
            # Compute next-period capital using bk_1d
            k[i+1]=np.matmul(x[i], bk_1d)
        
        x=np.asarray(x)
        # Compute time series of consumption
        c= np.multiply(np.power(k[0:T], alpha),a)+ (1-delta)*k[0:T]-k[1:T+1]
        # Monte Carlo realizations defined in condition (44) JMM (2011)
        y= np.multiply(np.divide(beta*np.power(c[1:T],-gam),np.power(c[0:T-1], -gam))*(1-delta+np.multiply(alpha*np.power(k[1:T],(alpha-1)),a[1:T])),k[1:T])
        # Compute a unit-free difference between the series from two iterations from JMM condition (10)
        dif_GSSA_1d = np.mean(abs(1-np.divide(k,k_old)))
        if checker == 1:
            print(dif_GSSA_1d) #for values checking
        
        # Compute new coefficients of the capital policy function using the OLS
        bk_hat_1d= np.matmul(inv(np.matmul(np.transpose(x[0:T-1,:]),x[0:T-1,:])),np.matmul(np.transpose(x[0:T-1,:]),y[0:T-1]))
        # Update the coefficients of the capital policy function using damping
        bk_1d = (kdamp*bk_hat_1d+ (1-kdamp)*bk_1d.reshape(1,3)).reshape(3,1)
        # Store the capital series
        k_old=np.copy(k)

    return y



def GSSA_poly(T, a, z, d, PF, zb, RM, penalty, normalize, dif_GSSA_D, kdamp, alpha, beta, delta, k, gam, y, k_old, a1, IM, n_nodes, weight_nodes, checker = 0):
    '''
    Stage 2 of 1 agent GSSA.
    ---------
    Parameters:
        T(int): Simulation points.
        a(numpy array): Productivity levels.
        z(numpy array): Data points on which the polynomial basis functions must be constructed.
        d(int): Number of polynomial.
        PF(binary): Polynomial family. 0=Ordinary, 1=Hermite.
        zb(numpy array): Matrix of means and standard deviations of state variables.
        RM(int): Regression method. 1=OLS, 2=LS-SVD, 3=LAD-PP, 4=LAD-DP, 5=RLS-Tikhonov, 6=RLS-TSVD, 7=RLAD-PP.
        penalty(int): Degree of regularization for a regularization methods, in particular RM= 5~7. For RM=6 needs to be positive, others negative.
        normalize(binary): Option of normalizing the data. 0=unnormalized data, 1=normalized data.   
        dif_GSSA_D: The initial difference between the series from two iterations in the convergence criterion.
        kdamp(float): Damping parameter. 
        alpha(float): Capital share in output.
        beta(float): Discount factor.
        delta(float): Depreciation rate.
        k(numpy array): Capital series.
        gam(float): Utility-function parameter.
        y(numpy array): Monte Carlo realizations from stage 1 GSSA.
        k_old(numpy array): Initialize the series of next-period capital.
        a1(numpy array): next-period productivity levels for each integration node using JMM condition (4).
        IM(int): Integration method for solution computation. 0=a one-node Monte Carlo method (default), 0=a one-node Monte Carlo method (default). 1~10 are Gauss-Hermite quadrature rules with 1~10 nodes, respectively.
        n_nodes(int): Number of integration nodes.
        weight_nodes(numpy array): weight_nodes are integration weights for Gauss-Hermite quadrature integration rule with IM nodes.
        checker(binary): A checker for the while loop, used to validate if the loop is running correctly or not. Default value is 0.

    --------
    Returns:
        bk_D (numpy array): The coefficients of the capital policy.
    '''
    # Construct the matrix of explanatory variables X on the series of state variables from the previously computed time-series solution
    X = Ord_Herm_Pol_1(z=z, D=d, PF=PF, zb=zb)
    # Compute the initial guess on the coefficients using the chosen regression method
    bk_D = Num_Stab_Approx(X[0:T-1,:],y[0:T-1,:],RM,penalty,normalize)

    # The main iterative cycle of GSSA

    # Checking convergence
    while dif_GSSA_D > 10**(-4)/10**(d)*kdamp:
        # Generate time series of capital
        for i in range(T):
            # The basis functions of a polynomial of degree D at time t
            X[i]= Ord_Herm_Pol_1(np.array([k[i],float(a[i])]).reshape((1,2)), D=d, PF=PF, zb=zb)
            # Compute next-period capital using bk_D
            k[i+1] = np.matmul(X[i], bk_D)
        
        # Compute time series of consumption
        c = np.multiply(np.power(k[0:T], alpha).reshape((T,1)), a) + (1-delta)*k[0:T].reshape(T,1)-k[1:T+1].reshape(T,1)

        # Preparing containers
        k1 = k[1:T+1]
        c1 = np.zeros((T,n_nodes))

        # The one-node Monte Carlo integration method approximates y with a realization of the integrand in the next period
        if IM == 0:
            y = np.multiply(np.divide(beta*np.power(c[1:T],-gam),np.power(c[0:T-1], -gam))*(1-delta+np.multiply(alpha*np.power(k[1:T],(alpha-1)).reshape(T-1,1),a[1:T])),k[1:T].reshape(T-1, 1))
        else:
            # Deterministic integration methods approximate y as a weighted average of the values of the integrand in the given nodes with the given weights
            for i in range(n_nodes):
                k2_prep = np.concatenate((k1.reshape(T,1), a1[:,i].reshape(T,1)), 1)
                # For each t, k2 is capital of period t+2
                k2= np.matmul(Ord_Herm_Pol_1(k2_prep,d,PF,zb), bk_D)
                # For each t, c1 is next-period consumption
                c1[:,i] = np.ravel((np.multiply(np.power(k1, alpha),a1[:,i])+ (1-delta)*k1).reshape(T,1) - k2)
            # Duplicate k1 n_nodes times to create a matrix with n_nodes identical columns
            k1_dupl = np.tile(k1, (n_nodes,1)).transpose()
            # Duplicate c n_nodes times to create a matrix with n_nodes identical columns
            c_dupl = np.tile(np.ravel(c), (n_nodes,1)).transpose()
            # Condition (8) in JMM (2011)
            y = np.matmul(np.multiply(np.divide(beta*np.power(c1,-gam),np.power(c_dupl, -gam))*(1-delta+np.multiply(alpha*np.power(k1_dupl,(alpha-1)),a1)),k1_dupl), weight_nodes)
        
        # The convergence criterion is adjusted to the damping parameter
        dif_GSSA_D = np.mean(abs(1-np.divide(k,k_old)))
        if checker == 1:
            print(dif_GSSA_D) #for values checking
        # Compute new coefficients of the capital policy function using the chosen approximation method
        bk_hat_D = Num_Stab_Approx(X[0:T-1,:],y[0:T-1,:],RM,penalty,normalize)
        # Update the coefficients of the capital policy function using damping
        bk_D = (kdamp*bk_hat_D+ (1-kdamp)*bk_D)
        # Store the capital series 
        k_old=np.copy(k)

    return bk_D


def GSSA_ShowcaseResult():
    '''
    This function aims to simply showcase the python implementation of GSSA with respect to the authors original Matlab codes.
    ---------
    Notice: All values are predetermined as the original matlab codes, also the comments closely follow the original
    comments to ensure transparency of the translation.

    ----------
    Output:
    showcase_result(Pandas Dataframe)
    '''
    ############################################
    #Roadmap of GSSA:
    #

    ############################################

    #We wont be simulating any data, we will use the data provided to ensure the result
    df =reading("epsi10000.csv")

    #1. Choose the simulation length
    #--------------------------------

    #Choose the simulation length for the solution procedure, T<=10,000
    T  = 10000

    #2. Model's parameters

    # Utility-function parameter                              
    gam     = 1   
    # Capital share in output     
    alpha   = 0.36 
    # Discount factor
    beta    = 0.99  
    # Depreciation rate   
    delta   = 0.02 
    # Persistence of the log of the productivity level     
    rho     = 0.95  
    # Standard deviation of shocks to the log of the productivity level   
    sigma   = 0.01 

    # 3. Steady state of capital (via Euler equation)
    ks = ( (1-beta+beta*delta) / (alpha*beta) )**(1/(alpha-1) )

    # 4. Initial condition
    k = np.array([ks]*(T+1))
    a= [1]*(T)
    epsi = df.iloc[:,0].astype(float)*sigma
    epsi=epsi.tolist()
    for i in range(1, T):
        a[i]=a[i-1]**(rho)*math.exp(epsi[i])
    a=np.asarray(a)
    kdamp = 0.01    
    dif_GSSA_1d = 1e+10  
    bk_1d  = np.array([0., 0.95, ks*0.05])
    bk_1d= np.reshape(bk_1d, (3,1))
    k_old = [ks+1]*(T+1)
    #GSSA stage 1 Initial guess by implementing the Monte Carlo simulation
    start = time.time()
    y= GSSA_main_cycle(T, gam, alpha, beta, delta, kdamp, dif_GSSA_1d, a, bk_1d, k_old, k)
    end = time.time()
    elapsed_time = end-start
    y = y.reshape((y.shape[0],1)) #make sure y is in the right shape
    #The GSSA parameters
    kdamp = 0.1
    dif_GSSA_D = 1e+10
    #The matrices of the polynomial coefficients
    D_max  = 5 #because of python
    npol = np.array([3, 6, 10, 15, 21])

    # 13. Choose an integration method for computing solutions  
    IM  = 10
    n_nodes,epsi_nodes, weight_nodes= GH_Quadrature(Qn=10, N=1, vcv=sigma**2)

    #make sure to change a into the right shape
    a = np.reshape(a, (T, 1))
    a1 = np.matmul(np.power(a,rho), np.exp(epsi_nodes.transpose()))

    #14. Choose a regression specification 
    RM = 6           # Choose a regression method: 
                 # 1=OLS,          2=LS-SVD,   3=LAD-PP,  4=LAD-DP, 
                 # 5=RLS-Tikhonov, 6=RLS-TSVD, 7=RLAD-PP, 8=RLAD-DP
    normalize = 1    # Option of normalizing the data; 0=unnormalized data; 
                 # 1=normalized data                    
    penalty = 7      # Degree of regularization for a regularization methods, 
                 # RM=5,6,7,8 (must be negative, e.g., -7 for RM=5,7,8 
                 # and must be positive, e.g., 7, for RM=6)
    PF = 0           # Choose a polynomial family; 0=Ordinary (default); # 1=Hermite
    # 15. Initialize the capital series
    zb = np.matrix([[np.mean(k[0:T]), np.mean(a[0:T])], [np.std(k[0:T]), np.std(a[0:T])]])
    z = np.concatenate((k[0:T].reshape(T,1), a[0:T].reshape(T,1)), axis=1)
    k_old = [ks+1]*(T+1)

    #
    BK = []
    Time = []
    for d in range(1, D_max+1):
        start = time.time()
        BK.append(GSSA_poly(T, a, z, d, PF, zb, RM, penalty, normalize, dif_GSSA_D, kdamp, alpha, beta, delta, k, gam, y, k_old, a1, IM, n_nodes, weight_nodes, checker= 0))
        end = time.time()
        Time.append(end-start)
    #Accuracy test:

    # we also will also just use the given data the authors generated
    # Notice for the real simulation we will differ from here
    T_test = 10200
    df =reading("epsi_test.csv")
    epsi_test = sigma*df.to_numpy().astype(float)
    a_test = [1]
    for i in range(1,T_test):
        value = a_test[i-1]**(rho)*math.exp(float(epsi_test[i]))
        a_test.append(value)

    IM_test = 10
    k_test = [ks]

    #construct containers for results
    result_max = []
    result_mean = []
    result_time = []
    for d in range(1, D_max+1):
        #refressing k_test to make sure that k_test is always 10200
        #k_test = [ks]
        for i in range(T_test):
            X_test = Ord_Herm_Pol_1(np.array([k_test[i], a_test[i]]).reshape([1,2]),d,PF,zb) # D = 1 for now, we will plug this in another for loop
            value = float(np.matmul(X_test, BK[d-1]))
            k_test.append(value)

        # testing it below
        discard = 200 #new defined value
        mean_error, max_error, error_time = Accuracy_Test_1(sigma,rho,beta,gam,alpha,delta,k_test,a_test,BK[d-1],d,IM_test,PF,zb,discard)
        result_max.append(max_error)
        result_mean.append(mean_error)
        result_time.append(error_time)

    #construct df
    showcase_result = pd.DataFrame(result_max, columns=["Maximum Error"])
    showcase_result["Mean Error"] = result_mean
    showcase_result["Time"] = Time
    showcase_result["Error Time"] = result_time
    showcase_result["Polynomial Degree"]= list(range(1,D_max+1))
    showcase_result["Total Time"] = showcase_result["Error Time"] + showcase_result["Time"]
    showcase_result["Rounded Total Time"]=showcase_result["Total Time"].round(decimals=2)
    showcase_result["Original Mean Error"]= 10**(showcase_result["Mean Error"])
    showcase_result["Original Max Error"]= 10**(showcase_result["Maximum Error"])
    showcase_result.set_index(["Polynomial Degree"])
    
    return showcase_result


def GSSA_1_agent(DC=5, IM=10, RM=6 ,normalize=1, penalty=7, PF=0):
    '''
    This is a modified version of the showcase code, notice that the total simulations are unified to 3000 following JMM (2011).
    --------
    Arguments:
        DC(int): Default is 5
        IM(int): Default is 10
        RM(int): Default is 6
        normalize(binary): Default is 1
        penalty(int): Default is 7
        PF(binary):  Default is 0

    -----------
    Output:
        result_max(list):
        result_mean(list):
        error_time(list):
        BK(list): 
        Time(list):
        stage1_time(float): 
    '''

    #################
    #GSSA stage 1
    #################

    # New simulation length
    # We will be using 3000 similar to the paper
    T = 3000
    #set seed to ensure reproduction
    np.random.seed(123)
    df_prep = randn2(T)
    df = pd.DataFrame(df_prep)

    # All parameters are the same as the original codes
    # Please consult the showcase code for explanations
    gam = 1       
    alpha = 0.36     
    beta = 0.99     
    delta = 0.02     
    rho = 0.95     
    sigma = 0.01    
    ks = ( (1-beta+beta*delta) / (alpha*beta) )**(1/(alpha-1) )
    k = np.array([ks]*(T+1))
    a= [1]*(T)
    epsi = df.iloc[:,0].astype(float)*sigma
    epsi=epsi.tolist()
    for i in range(1, T):
        a[i]=a[i-1]**(rho)*math.exp(epsi[i])
    a=np.asarray(a)
    kdamp = 0.01    
    dif_GSSA_1d = 1e+10  
    bk_1d  = np.array([0., 0.95, ks*0.05])
    bk_1d= np.reshape(bk_1d, (3,1))
    k_old = [ks+1]*(T+1)
    start = time.time()
    y= GSSA_main_cycle(T, gam, alpha, beta, delta, kdamp, dif_GSSA_1d, a, bk_1d, k_old, k)
    end = time.time()
    stage1_time = end-start

    #make sure y is in the right shape
    y = y.reshape((y.shape[0],1))

    ####################
    #GSSA stage 2
    ####################
    #The GSSA parameters
    kdamp = 0.1
    dif_GSSA_D = 1e+10
    #The matrices of the polynomial coefficients
    D_max  = DC 
    npol = np.array([3, 6, 10, 15, 21])

    #Choosing integration method, please consult showcase for the specifics 
    IM  = IM
    n_nodes,epsi_nodes, weight_nodes= GH_Quadrature(Qn=10, N=1, vcv=sigma**2)
    a = np.reshape(a, (T, 1))
    a1 = np.matmul(np.power(a,rho), np.exp(epsi_nodes.transpose()))

    #Please consult the showcase code for the existing methods
    RM = RM               
    normalize = normalize                      
    penalty = pt      
    PF = PF          
    zb = np.matrix([[np.mean(k[0:T]), np.mean(a[0:T])], [np.std(k[0:T]), np.std(a[0:T])]])
    z = np.concatenate((k[0:T].reshape(T,1), a[0:T].reshape(T,1)), axis=1)
    k_old = [ks+1]*(T+1)
    BK = []
    Time = []
    for d in range(1, D_max+1):
        start = time.time()
        BK.append(GSSA_poly(T, a, z, d, PF, zb, RM, penalty, normalize, dif_GSSA_D, kdamp, alpha, beta, delta, k, gam, y, k_old, a1, IM, n_nodes, weight_nodes, checker= 0))
        end = time.time()
        Time.append(end-start)

    ########################
    #Accuracy testing
    ########################
    T_test = 10000
    np.random.seed(100)
    df_prep = randn2(T_test)
    df = pd.DataFrame(df_prep)
    epsi_test = sigma*df.to_numpy().astype(float)
    a_test = [1]
    for i in range(1,T_test):
        value = a_test[i-1]**(rho)*math.exp(float(epsi_test[i]))
        a_test.append(value)

    IM_test = IM
    k_test = [ks]
    result_max = []
    result_mean = []
    result_time = []
    for d in range(1, D_max+1):
        for i in range(T_test):
            X_test = Ord_Herm_Pol_1(np.array([k_test[i], a_test[i]]).reshape([1,2]),d,PF,zb)
            value = float(np.matmul(X_test, BK[d-1]))
            k_test.append(value)
        
        #The T_test - discard needs to match the exact numbers of simulations
        discard = 7000 
        mean_error, max_error, error_time = Accuracy_Test_1(sigma,rho,beta,gam,alpha,delta,k_test,a_test,BK[d-1],d,IM_test,PF,zb,discard)
        result_max.append(max_error)
        result_mean.append(mean_error)
        result_time.append(error_time)
    
    return result_max, result_mean, error_time, BK, Time, stage1_time