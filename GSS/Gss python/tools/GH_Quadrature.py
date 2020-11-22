import numpy as np
import math
from numpy.linalg import cholesky

#This is the python implementation equivalent of GH_Quadrature by Lilia Maliar and Serguei Maliar, this version is constructed by Viktor Cheng which included multiple optimisations under the python framework.

def GH_Quadrature(N,vcv,Qn=10):
    ''''''
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
    else: #The default option is given Qn=10
        eps = np.array([3.436159118837738,2.532731674232790,1.756683649299882,1.036610829789514,0.3429013272237046,-0.3429013272237046,-1.036610829789514,-1.756683649299882,-2.532731674232790,-3.436159118837738])
        weight = np.array([7.640432855232621e-06,0.001343645746781233,0.03387439445548106,0.2401386110823147,0.6108626337353258,0.6108626337353258,0.2401386110823147,0.03387439445548106,0.001343645746781233,7.640432855232621e-06])
    
    n_nodes = Qn**N
    z1 = np.ones([n_nodes,N]).astype(float)
    w1i = np.ones([n_nodes,N]).astype(float)

    for i in range(N):
        z1[:,i]=np.tile(np.repeat(eps,10**(i)), 10**(N-i-1))
        w1i[:,i]=np.tile(np.repeat(weight,10**(i)), 10**(N-i-1))
    w1=w1i[:,0]*w1i[:,1]*w1i[:,2]
    z= math.sqrt(2)*z1
    weight_nodes= w1/(math.sqrt(math.pi)**N)
    epsi_nodes = np.matmul(z,cholesky(vcv))

    return n_nodes,epsi_nodes,weight_nodes


