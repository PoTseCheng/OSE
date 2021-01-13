import time
import numpy as np
from tools.Ord_Herm_Pol_1 import*
from tools.GH_Quadrature import*


def Accuracy_Test_1(sigma,rho,beta,gam,alpha,delta,k,a,bk,D,IM,PF,zb,discard):
    '''
    '''
    start = time.time()
    n_nodes,epsi_nodes, weight_nodes = GH_Quadrature(Qn=IM, N=1, vcv=sigma**2)
    #initialise empty list
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