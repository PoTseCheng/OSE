import numpy as np
from tools.GH_Quadrature import*
from tools.Ord_Herm_Pol_1 import*
from tools.Num_Stab_Approx import*




def GSSA_main_cycle(T, gam, alpha, beta, delta, kdamp, dif_GSSA_1d, a, bk_1d, k_old, k, checker = 0):
    '''
    Compute a first-degree polynomial solution using the one-node Monte Carlo integration method.

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
    Compute the polynomial solution.

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


    #to do: Maybe make the function easier, in the sense allow less inputs