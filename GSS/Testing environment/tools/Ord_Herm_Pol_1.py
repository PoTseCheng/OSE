import numpy as np




def Ord_Herm_Pol_1(z,D,PF,zb):
    '''
    
    '''
    # assumed the imported matrix is a numpy matrix### like above
    # construct container for results
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