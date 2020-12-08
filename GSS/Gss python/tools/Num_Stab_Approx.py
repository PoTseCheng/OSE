import pandas as pd
import numpy as np
from numpy.linalg import inv
from numpy.linalg import svd
from scipy.optimize import linprog


def Num_Stab_Approx(x,y,RM,penalty,normalised):
    '''
    
    '''
    T = x.shape[0]
    n = x.shape[1]
    N = y.shape[1] #Compute the number of columns in Y
    
    if normalised == 1 or RM >= 5:
        X1 = np.divide(
            (x[:,1:n] - np.matmul(np.ones((T, 1)), x[:,1:n].mean(axis=0).reshape(1, n-1))),
            np.matmul(np.ones((T, 1)), np.array([np.std(x[:,i], ddof=1) for i in range(1,n)]).reshape(1,n-1))
            )#unfortunately np.std seems cannot return an array for each column, therefore I was forced to use a list comprehension here
        
        X1 = X1.astype(float)
        Y1 = np.divide(
            (y - np.ones((T,1))*np.mean(y)),
            np.matmul(np.ones((T,1)), np.std(y, ddof=1).reshape(1,1))
            )
        Y1 = Y1.astype(float)
        n1 = n-1 # Number of coefficients in a regression with normalized data is reduced by 1 (no intercept)

        # regression methods
        if RM == 1:
            B = np.matmul(
                inv(np.matmul(X1.transpose(),X1).astype(float)), #need to convert the output to float or else return error
                np.matmul(X1.transpose(),Y1))
        # LS-SVD
        elif RM == 2:
            U, S, Vh = svd(X1, full_matrices=False)
            V = Vh.T
            S_inv = np.diag(1/S)
            B = np.matmul(np.matmul(np.matmul(V,S_inv),U.transpose()),Y1) 
        elif RM == 3:
            # This method is strongly discouraged! It is one of the slowest among all other regression construct the boundaries
            BND = [(-100, 100)]*n1 + [(0, float("inf"))]*2*T
            f = [0]*n1+ [1]*2*T
            #specify the Aeq an beq
            Aeq = np.concatenate((X1, np.identity(T), -np.identity(T)), axis=1)
            B =[]
            #solve the equation
            for j in range(N):
                beq = Y1[:,j].tolist()
                result = linprog(f, A_eq = Aeq, b_eq = beq, bounds= BND , method='revised simplex')
                B.append(result.x) 
                # Or so it should. Unfortunately, the Scipy package will took ages to even finish computing one optimisation solution given the size of Aeq and beq. Without other choices, I import the function that does this part of computation from matlab, which is able to finish the computation within seconds.
        
        elif RM == 4:
            # Also extremely slow
            # Define Boundary
            BND = [(-1,1)]*T
            # Define aeq and beq for linprog
            Aeq = X1.transpose()
            beq = np.zeros((n1, 1))
            B = []
            for i in range(N):
                f = -Y1[:,i]
                result = linprog(f, A_eq = Aeq, b_eq = beq, bounds= BND)
                B.append(result.x)
        
        # RLS-Tikhonov
        elif RM == 5:
            B = np.matmul(
                inv(np.matmul(X1.transpose(),X1) + (T/n1)*np.identity(n1)*10**(penalty)),
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

        # normalised the output
        B2 = np.multiply((1/np.asarray([np.std(x[:,i], ddof=1) for i in range(1,n)])).reshape((2,1))*np.std(y, ddof=1),B)
        B1 = (y.mean()- np.matmul(np.asarray([x[:,i].mean() for i in range(1,n)]),B2)).reshape(1,1)
        B = np.concatenate((B1,B2))
        return B
        
        
    else:
        X1 = x
        Y1 = Y          # Leave Y without changes
        n1 = n
        # Regression methods
        
        if RM == 1:
            B = np.matmul(
                inv(np.matmul(X1.transpose(),X1).astype(float)), #need to convert the output to float or else return error
                np.matmul(X1.transpose(),Y1)
                )
            return B
        
        #LS-SVD
        elif RM == 2:
            U, S, Vh = svd(X1, full_matrices=False)
            V = Vh.T
            S_inv = np.diag(1/S)
            B = np.matmul(np.matmul(np.matmul(V,S_inv),U.transpose()),Y1)
            return B
        
        elif RM == 3:
            # #This method is strongly discouraged! It is one of the slowest among all other regression
            # construct the boundaries
            BND = [(-100, 100)]*n1 + [(0, float("inf"))]*2*T
            f = [0]*n1+ [1]*2*T
            # specify the Aeq an beq
            Aeq = np.concatenate((X1, np.identity(T), -np.identity(T)), axis=1)
            B =[]
            #solve the equation
            
            for j in range(N):
                beq = Y1[:,j].tolist()
                result = linprog(f, A_eq = Aeq, b_eq = beq, bounds= BND)
                B.append(result.x) 
            #Or so it should. Unfortunately, the Scipy package will took ages to even finish computing one optimisation solution given the size of Aeq and beq. Without other choices, I import the function that does this part of computation from matlab, which is able to finish the computation within seconds.
            return B
        
        elif RM == 4:
            # Also extremely slow
            # Define Boundary
            BND = [(-1,1)]*T
            #Define aeq and beq for linprog
            Aeq = X1.transpose()
            beq = np.zeros((n1, 1))
            B = []
            for i in range(N):
                f = -Y1[:,i]
                result = linprog(f, A_eq = Aeq, b_eq = beq, bounds= BND)
                B.append(result.x)
            return B
                
    