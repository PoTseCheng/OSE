import pandas as pd
import numpy as np
from numpy.linalg import inv
from numpy.linalg import svd
from scipy.optimize import linprog

def Num_Stab_Approx(X,Y,RM,penalty,normalised):
    '''

    '''
    T = x.shape[0]
    n = x.shape[1]
    N = y.shape[1] #Compute the number of columns in Y

    if normalised == 1 or RM >= 5:
        X1 = np.divide(
            (x[:,1:n] - np.matmul(np.ones((T, 1)), x[:,1:n].mean(axis=0).reshape(1, n-1))),
            np.matmul(np.ones((T, 1)), np.array([np.std(x[:,i], ddof=1) for i in range(1,n)]).reshape(1,n-1))
        )
        X1 = X1.astype(float)
    #unfortunately np.std seems cannot return an array for each column, therefore I was forced to use a list comprehension here
    
        Y1 = np.divide(
            (y - np.ones((T,1))*np.mean(y)),
            np.matmul(np.ones((T,1)), np.std(y, ddof=1).reshape(1,1))
            )
        Y1 = Y1.astype(float)
        n1 = n-1 # Number of coefficients in a regression with normalized data is reduced by 1 (no intercept)
    else:
        X1 = x
        Y1 = y          # Leave Y without changes
        n1 = n
    
    # 3. Regression methods
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
    
    elif RM == 3:
        BND = [(-100, 100)]*n1 + [(0, None)]*2*T
        f = np.asarray([0]*n1+ [1]*2*T)
        Aeq = np.concatenate((X1, np.eye(T), -np.eye(T)), axis=1)
        B =[]
        #solve the equation
        for i in range(N):
            beq = Y1[:,i]
            result = linprog(f, A_eq = Aeq, b_eq = beq, bounds= BND, method="highs-ipm")
            B.append(list(result.x[0:n1]))
            B = np.asarray(B).reshape(n1,1)
    
    elif RM == 4:
        #This method is unavaliable
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
        Sr_inv = np.diag(np.divide(1., S))
        B = np.matmul(
                    np.matmul(
                        np.matmul(V, Sr_inv),
                        U.transpose()
                             ),
                    Y1)
    #LADPP
    elif RM == 7:
        #we can just use the default setting from scipy as the lower and upper will be the same
        f= np.r_[T/n1*10**(penalty)*np.ones(4,),np.ones(2*T,)]
        Aeq= np.c_[X1,-X1, np.eye(T), -np.eye(T)]
        B =[]
        #solve the equation
        for i in range(N):
            beq = Y1[:,i]
            result = linprog(f, A_eq = Aeq, b_eq = beq, method="highs-ipm")
            B.append(list(result.x[0:n1]-result.x[n1:2*n1]))
            B = np.asarray(B).reshape(n1,1)

    #RM == 8 is unavaliable(see notebook)
    if normalised == 1 or RM >= 5:
        B2 = np.multiply((1/np.asarray([np.std(x[:,i], ddof=1) for i in range(1,n)])).reshape((2,1))*np.std(y, ddof=1),B)
        B1 = (y.mean()- np.matmul(np.asarray([x[:,i].mean() for i in range(1,n)]),B2)).reshape(1,1)
        B = np.concatenate((B1,B2))
    else:
        B = np.matmul(
            inv(np.matmul(X1.transpose(),X1).astype(float)), #need to convert the output to float or else return error
            np.matmul(X1.transpose(),Y1)
                     )
    
    return B