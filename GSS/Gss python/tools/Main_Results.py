import numpy as np
from tools.GH_Quadrature import*
from tools.Ord_Herm_Pol_1 import*
from tools.Num_Stab_Approx_new import*
from tools.Accuracy_Test_1 import*
from tools.auxiliary import*
import math
import time




def GSSA_main_cycle(T, gam, alpha, beta, delta, kdamp, dif_GSSA_1d, a, bk_1d, k_old, k, checker = 0):
    '''
    Stage 1 of GSSA

            Parameters:
                    T (int): A decimal integer
                    gam (int): Another decimal integer
                    alpha:
                    beta:
                    delta: 
                    kdamp: 
                    dif_GSSA_1d: 
                    a: 
                    bk_1d:
                    k_old:
                    k:
                    checker (int): A checker for the while loop, used to validate if the loop is running correctly or not.

            Returns:
                    numpy array: A numpy array of floats with the shape of (9999,)
    '''
    
    while dif_GSSA_1d > 1e-4*kdamp:
        x=[]
        for i in range(0,T):
            x.append([1,k[i],a[i]])
            k[i+1]=np.matmul(x[i], bk_1d)
        x=np.asarray(x)
        c= np.multiply(np.power(k[0:T], alpha),a)+ (1-delta)*k[0:T]-k[1:T+1]
        y= np.multiply(np.divide(beta*np.power(c[1:T],-gam),np.power(c[0:T-1], -gam))*(1-delta+np.multiply(alpha*np.power(k[1:T],(alpha-1)),a[1:T])),k[1:T])
        dif_GSSA_1d = np.mean(abs(1-np.divide(k,k_old)))
        if checker == 1:
            print(dif_GSSA_1d) #for values checking
        bk_hat_1d= np.matmul(inv(np.matmul(np.transpose(x[0:T-1,:]),x[0:T-1,:])),np.matmul(np.transpose(x[0:T-1,:]),y[0:T-1]))
        bk_1d = (kdamp*bk_hat_1d+ (1-kdamp)*bk_1d.reshape(1,3)).reshape(3,1)
        k_old=np.copy(k)

    return y



def GSSA_poly(T, a, z, d, PF, zb, RM, penalty, normalize, dif_GSSA_D, kdamp, alpha, beta, delta, k, gam, y, k_old, a1, IM, n_nodes, weight_nodes, checker = 0):
    '''
    Stage 2 of GSSA in Judd et. al (2011) including accuracy test.

            Parameters:
                    T:
                    a:
                    z:
                    d:
                    PF:
                    zb:
                    RM:
                    penalty: 
                    normalize: 
                    dif_GSSA_D:
                    kdamp:
                    alpha:
                    delta:
                    k:
                    gam:
                    y:
                    k_old:
                    a1:
                    IM:
                    n_nodes:
                    weight_nodes:
                    bk_D:
                    checker (int): A checker for the while loop, used to validate if the loop is running correctly or not. Default value is 0.

            Returns:
                    bk_D (numpy array): 
    '''
    X = Ord_Herm_Pol_1(z=z, D=d, PF=PF, zb=zb)
    bk_D = Num_Stab_Approx(X[0:T-1,:],y[0:T-1,:],RM,penalty,normalize)
    while dif_GSSA_D > 10**(-4)/10**(d)*kdamp:
        for i in range(T):
            X[i]= Ord_Herm_Pol_1(np.array([k[i],float(a[i])]).reshape((1,2)), D=d, PF=PF, zb=zb)
            k[i+1] = np.matmul(X[i], bk_D)
        c = np.multiply(np.power(k[0:T], alpha).reshape((T,1)), a) + (1-delta)*k[0:T].reshape(T,1)-k[1:T+1].reshape(T,1)
        k1 = k[1:T+1]
        c1 = np.zeros((T,n_nodes))
        if IM == 0:
            y = np.multiply(np.divide(beta*np.power(c[1:T],-gam),np.power(c[0:T-1], -gam))*(1-delta+np.multiply(alpha*np.power(k[1:T],(alpha-1)).reshape(T-1,1),a[1:T])),k[1:T].reshape(T-1, 1))
        else:
            for i in range(n_nodes):
                k2_prep = np.concatenate((k1.reshape(T,1), a1[:,i].reshape(T,1)), 1)  #preparation for k2 
                k2= np.matmul(Ord_Herm_Pol_1(k2_prep,d,PF,zb), bk_D)  #1 equals D
                c1[:,i] = np.ravel((np.multiply(np.power(k1, alpha),a1[:,i])+ (1-delta)*k1).reshape(T,1) - k2)
            k1_dupl = np.tile(k1, (n_nodes,1)).transpose()
            c_dupl = np.tile(np.ravel(c), (n_nodes,1)).transpose()
            y = np.matmul(np.multiply(np.divide(beta*np.power(c1,-gam),np.power(c_dupl, -gam))*(1-delta+np.multiply(alpha*np.power(k1_dupl,(alpha-1)),a1)),k1_dupl), weight_nodes)
        dif_GSSA_D = np.mean(abs(1-np.divide(k,k_old)))
        if checker == 1:
            print(dif_GSSA_D) #for values checking
        bk_hat_D = Num_Stab_Approx(X[0:T-1,:],y[0:T-1,:],RM,penalty,normalize)
        bk_D = (kdamp*bk_hat_D+ (1-kdamp)*bk_D)
        k_old=np.copy(k)

    return bk_D


def GSSA_MainResult(D_Max=5, Rm=6, Norml=1, Pt=7, Pf=0 ):
    '''
    This function aims to showcase the important features of GSSA, namely comparing different regression methods, while holding 
    other arguments fixed. This is to avoid a super long input as GSSA provides huge degree of freedom for defining the model.

    Arguments:
    D_Max(int): Maximum degree of a polynomial, can be 1 to 5.
    Rm(int):Regression methods:
    Norml(int): Normalised the data or not, by default is 1. 
    Pt(int): Degree of regularization, needs to match with certain Rm,
    Pf(int):

    Output:
    Pandas Dataframe

    
    '''
    #todo:
    #add more comments

    ############################################
    #Roadmap of GSSA:
    #

    ############################################
    df =reading("epsi10000.csv")
    T  = 10000           #Choose the simulation length for the solution procedure, T<=10,000                    
    gam     = 1        # Utility-function parameter
    alpha   = 0.36     # Capital share in output
    beta    = 0.99     # Discount factor
    delta   = 0.02     # Depreciation rate 
    rho     = 0.95     # Persistence of the log of the productivity level
    sigma   = 0.01    # Standard deviation of shocks to the log of the productivity level
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
    #Main GSSA cycle
    start = time.time()
    y= GSSA_main_cycle(T, gam, alpha, beta, delta, kdamp, dif_GSSA_1d, a, bk_1d, k_old, k)
    end = time.time()
    elapsed_time = end-start
    y = y.reshape((y.shape[0],1)) #make sure y is in the right shape
    #The GSSA parameters
    kdamp = 0.1
    dif_GSSA_D = 1e+10
    #The matrices of the polynomial coefficients
    D_max  = D_Max #because of python
    npol = np.array([3, 6, 10, 15, 21])

    # 13. Choose an integration method for computing solutions  
    IM  = 10
    n_nodes,epsi_nodes, weight_nodes= GH_Quadrature(Qn=10, N=1, vcv=sigma**2)

#make sure to change a into the right shape
    a = np.reshape(a, (T, 1))
    a1 = np.matmul(np.power(a,rho), np.exp(epsi_nodes.transpose()))

    #14. Choose a regression specification 
    RM = Rm           # Choose a regression method: 
                 # 1=OLS,          2=LS-SVD,   3=LAD-PP,  4=LAD-DP, 
                 # 5=RLS-Tikhonov, 6=RLS-TSVD, 7=RLAD-PP, 8=RLAD-DP
    normalize = Norml    # Option of normalizing the data; 0=unnormalized data; 
                 # 1=normalized data                    
    penalty = Pt      # Degree of regularization for a regularization methods, 
                 # RM=5,6,7,8 (must be negative, e.g., -7 for RM=5,7,8 
                 # and must be positive, e.g., 7, for RM=6)
    PF = Pf           # Choose a polynomial family; 0=Ordinary (default); # 1=Hermite
    # 15. Initialize the capital series
    zb = np.matrix([[np.mean(k[0:T]), np.mean(a[0:T])], [np.std(k[0:T]), np.std(a[0:T])]])
    z = np.concatenate((k[0:T].reshape(T,1), a[0:T].reshape(T,1)), axis=1)
    k_old = [ks+1]*(T+1)

    #Solution:
    BK = []
    Time = []
    for d in range(1, D_max+1):
        start = time.time()
        BK.append(GSSA_poly(T, a, z, d, PF, zb, RM, penalty, normalize, dif_GSSA_D, kdamp, alpha, beta, delta, k, gam, y, k_old, a1, IM, n_nodes, weight_nodes, checker= 0))
        end = time.time()
        Time.append(end-start)
    #Accuracy test:
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
