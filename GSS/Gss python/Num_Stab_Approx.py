# Generated with SMOP  0.41
from libsmop import *
# .\Num_Stab_Approx.m

    # Num_Stab_Approx.m is a routine that implements the approximation methods 
# described in "Numerically Stable and Accurate Stochastic Simulation 
# Approaches for Solving Dynamic Economic Models" by Kenneth L. Judd, Lilia 
# Maliar and Serguei Maliar, (2011), Quantitative Economics 2/2, 173–210 
# (henceforth, JMM, 2011).
    
    # This version: July 14, 2011. First version: August 27, 2009.
# -------------------------------------------------------------------------
# Inputs:  "X" is a matrix of dependent variables in a regression, T-by-n,
#          where n corresponds to the total number of coefficients in the 
#          original regression (i.e. with unnormalized data);
#          "Y" is a matrix of independent variables, T-by-N; 
#          "RM" is the regression (approximation) method, RM=1,...,8:  
#          1=OLS,          2=LS-SVD,    3=LAD-PP,    4=LAD-DP, 
#          5=RLS-Tikhonov, 6=RLS-TSVD,  7=RLAD-PP,   8=RLAD-DP;
#          "penalty"  is a parameter determining the value of the regulari-
#          zation parameter for a regularization methods, RM=5,6,7,8;
#          "normalize"  is the option of normalizing the data, 
#          0=unnormalized data,  1=normalized data
    
    # Outputs: "B" is a matrix of the regression coefficients 
# -------------------------------------------------------------------------
# Copyright © 2011 by Lilia Maliar and Serguei Maliar. All rights reserved. 
# The code may be used, modified and redistributed under the terms provided 
# in the file "License_Agreement.txt".
# -------------------------------------------------------------------------
    
    
@function
def Num_Stab_Approx(X=None,Y=None,RM=None,penalty=None,normalize=None,*args,**kwargs):
    varargin = Num_Stab_Approx.varargin
    nargin = Num_Stab_Approx.nargin

    # 1. Compute the dimensionality of the data
#------------------------------------------
    T,n=size(X,nargout=2)
# .\Num_Stab_Approx.m:32
    
    # number of regression coefficient to be computed
    N=size(Y,2)
# .\Num_Stab_Approx.m:34
    
    # to the total number of regressions to be ran
    
    # 2. Normalize the data
#----------------------
    if (normalize == 1) or (RM >= 5):
        # regularization method, ...
        X1=(X(arange(),arange(2,n)) - dot(ones(T,1),mean(X(arange(),arange(2,n))))) / (dot(ones(T,1),std(X(arange(),arange(2,n)))))
# .\Num_Stab_Approx.m:42
        Y1=(Y - dot(ones(T,1),mean(Y))) / (dot(ones(T,1),std(Y)))
# .\Num_Stab_Approx.m:44
        n1=n - 1
# .\Num_Stab_Approx.m:45
        # normalized data is reduced by 1 (no intercept)
    else:
        X1=copy(X)
# .\Num_Stab_Approx.m:48
        Y1=copy(Y)
# .\Num_Stab_Approx.m:49
        n1=copy(n)
# .\Num_Stab_Approx.m:50
    
    # 3. Regression methods
#----------------------
    
    # 3.1 OLS
#--------
    if RM == 1:
        B=dot(dot(inv(dot(X1.T,X1)),X1.T),Y1)
# .\Num_Stab_Approx.m:59
        # estimator; note that all the regressions,
                         # j=1,...,N, are ran at once
        # 3.2 LS-SVD
#-----------
    else:
        if RM == 2:
            U,S,V=svd(X1,0,nargout=3)
# .\Num_Stab_Approx.m:67
            # "economy size" SVD in MATLAB; matrices U, V and 
                         # S are defined by X1=U*S*V', where U is T-by-n1, 
                         # S is n1-by-n1, and V is n1-by-n1
            S_inv=diag(1.0 / diag(S(arange(1,n1),arange(1,n1))))
# .\Num_Stab_Approx.m:71
            B=dot(dot(dot(V,S_inv),U.T),Y1)
# .\Num_Stab_Approx.m:73
            # estimator (20) in JMM (2011)
            # For each regression j, RM 3, 4, 7, 8, we  solve a linear programming 
 # problem written in MATLAB:
 #                       min  f'*xlp 
 #                       s.t. Aeq*xlp=beq       (*)
 #                            Aineq*xlp<=bineq  (**) 
 #                            LB=xlp<=UB
 # where xlp is the vector of unknowns; Aeq is the matrix of coefficients 
 # in the equality restriction (*); Aineq is the matrix of coefficients in the
 # inequality restriction (**); beq and bineq are vectors of constants in the 
 # right sides of (*) and (**), respectively; LB means "lower bound", UB  
 # means "upper bound"
            # 3.3 LAD-PP
 #-----------
        else:
            if RM == 3:
                # We solve the linear programming problem (27)-(29) in JMM (2011). 
   # xlp=[B(:,j); ups_plus; ups_minus] where B(:,j) is the vector of coeffi-
   # cients in the regression j, ups_plus and ups_minus are the deviations; 
   # f'=[0,...,0,1,...,1] where zeros correspond to the coefficients in the 
   # objective function on B and ones correspond to the coefficients on 
   # ups_plus and ups_minus
                LB=concat([[zeros(n1,1) - 100],[zeros(dot(2,T),1)]])
# .\Num_Stab_Approx.m:97
                # -100, and lower bounds on ups_plus and ups_minus 
                         # are set to zero
                UB=concat([[zeros(n1,1) + 100],[inf(dot(2,T),1)]])
# .\Num_Stab_Approx.m:101
                # 100, and upper bounds on ups_plus and ups_minus 
                         # are set to infinity
                f=concat([[zeros(n1,1)],[ones(dot(2,T),1)]])
# .\Num_Stab_Approx.m:105
                Aeq=concat([X1,eye(T,T),- eye(T,T)])
# .\Num_Stab_Approx.m:106
                B=zeros(size(X1,2),N)
# .\Num_Stab_Approx.m:107
                # contain coefficients of all regressions; n1-by-N
                for j in arange(1,N).reshape(-1):
                    beq=Y1(arange(),j)
# .\Num_Stab_Approx.m:110
                    xlp,fval,exitflag,output,lambda_=linprog(f,[],[],Aeq,beq,LB,UB,[],nargout=5)
# .\Num_Stab_Approx.m:111
                    B[arange(),j]=xlp(arange(1,n1),1)
# .\Num_Stab_Approx.m:113
                    # regression j, xlp(1:n1,1), into the matrix B
                # 3.4 LAD-DP
 #-----------
            else:
                if RM == 4:
                    # We solve the linear programming problem (30)-(32) in JMM (2011). 
 # xlp=[q] where q is a vector of unknowns in (30)-(32) of JMM (2011)
                    LB=- ones(1,T)
# .\Num_Stab_Approx.m:122
                    UB=ones(1,T)
# .\Num_Stab_Approx.m:123
                    Aeq=X1.T
# .\Num_Stab_Approx.m:124
                    beq=zeros(n1,1)
# .\Num_Stab_Approx.m:125
                    B=zeros(size(X1,2),N)
# .\Num_Stab_Approx.m:126
                    # n1-by-N
                    for j in arange(1,N).reshape(-1):
                        f=- Y1(arange(),j)
# .\Num_Stab_Approx.m:130
                        # appears because (30)-(32) in JMM (2011) is a 
                         # minimization problem)
                        xlp,fval,exitflag,output,lambda_=linprog(f,[],[],Aeq,beq,LB,UB,[],nargout=5)
# .\Num_Stab_Approx.m:133
                        # lambda, on all the constraints
                        B[arange(),j]=lambda_.eqlin
# .\Num_Stab_Approx.m:136
                        # are equal to the Lagrange multipliers on the  
                         # equality constraint (*)
                    # 3.5 RLS-Tikhonov
#-----------------
                else:
                    if RM == 5:
                        B=dot(dot(inv(dot(X1.T,X1) + dot(dot(T / n1,eye(n1)),10 ** penalty)),X1.T),Y1)
# .\Num_Stab_Approx.m:145
                        # regularization parameter is T/n1*10^-penalty
                        # 3.6 RLS-TSVD
#-------------
                    else:
                        if RM == 6:
                            U,S,V=svd(X1,0,nargout=3)
# .\Num_Stab_Approx.m:152
                            # which is "economy size" SVD in MATLAB; matrices 
                         # U, V and S are defined by X1=U*S*V', where U is 
                         # T-by-n1, S is n1-by-n1, and V is n1-by-n1
                            r=sum((max(diag(S)) / diag(S)) <= 10 ** penalty)
# .\Num_Stab_Approx.m:157
                            # <= than a threshold level equal to 10^penalty
                            Sr_inv=zeros(n1)
# .\Num_Stab_Approx.m:160
                            Sr_inv[arange(1,r),arange(1,r)]=diag(1.0 / diag(S(arange(1,r),arange(1,r))))
# .\Num_Stab_Approx.m:160
                            B=dot(dot(dot(V,Sr_inv),U.T),Y1)
# .\Num_Stab_Approx.m:162
                            # the RLS-TSVD estimator (43) in JMM (2011)
                            # 3.7 RLAD-PP
#------------
                        else:
                            if RM == 7:
                                # We solve the linear programming problem (34)-(37) in JMM (2011). 
# xlp=[phi_plus; phi_minus; ups_plus; ups_minus] where phi_plus, phi_minus, ups_plus, ups_minus are defined in 
# (34)-(37) of JMM (2011)
                                LB=concat([[zeros(dot(2,n1),1)],[zeros(dot(2,T),1)]])
# .\Num_Stab_Approx.m:172
                                UB=[]
# .\Num_Stab_Approx.m:174
                                f=concat([[dot(dot(10 ** penalty,ones(dot(n1,2),1)),T) / n1],[ones(dot(2,T),1)]])
# .\Num_Stab_Approx.m:175
                                # (2*n1+2T)-by-1
                                Aeq=concat([X1,- X1,eye(T,T),- eye(T,T)])
# .\Num_Stab_Approx.m:178
                                B=zeros(size(X1,2),N)
# .\Num_Stab_Approx.m:180
                                # n1-by-N
                                for j in arange(1,N).reshape(-1):
                                    beq=Y1(arange(),j)
# .\Num_Stab_Approx.m:184
                                    xlp,fval,exitflag,output,lambda_=linprog(f,[],[],Aeq,beq,LB,UB,[],nargout=5)
# .\Num_Stab_Approx.m:185
                                    # and xlp(n1+1:2*n1,1) corresponds to phi_minus
                                    B[arange(),j]=xlp(arange(1,n1),1) - xlp(arange(n1 + 1,dot(2,n1)),1)
# .\Num_Stab_Approx.m:188
                                    # are given by the difference phi_plus - phi_minus;   
                         # store these coefficients into the matrix B
                                # 3.8 RLAD-DP
#------------
                            else:
                                if RM == 8:
                                    # We solve the linear programming problem (38)-(41) in JMM (2011). 
# xlp=[q] where q is a vector of unknowns in (38)-(41) of JMM (2011)
                                    LB=- ones(1,T)
# .\Num_Stab_Approx.m:199
                                    UB=ones(1,T)
# .\Num_Stab_Approx.m:200
                                    Aineq=concat([[X1.T],[- X1.T]])
# .\Num_Stab_Approx.m:201
                                    bineq=dot(dot(10 ** penalty,ones(dot(n1,2),1)),T) / n1
# .\Num_Stab_Approx.m:202
                                    # 2*n1-by-1
                                    B=zeros(size(X1,2),N)
# .\Num_Stab_Approx.m:205
                                    # n1-by-N
                                    for j in arange(1,N).reshape(-1):
                                        f=- Y1(arange(),j)
# .\Num_Stab_Approx.m:209
                                        # appears because (38)-(41) in JMM (2011) is a 
                         # minimization problem)
                                        xlp,fval,exitflag,output,lambda_=linprog(f,Aineq,bineq,[],[],LB,UB,[],nargout=5)
# .\Num_Stab_Approx.m:212
                                        # lambda, on all the constraints; phi_plus and phi_minus 
                         # (defined in (38)-(41) in JMM, 2011) are equal to  
                         # the Lagrange multipliers lambda.ineqlin(1:n1) 
                         # and lambda.ineqlin(n1+1:2*n1), respectively
                                        B[arange(),j]=lambda_.ineqlin(arange(1,n1)) - lambda_.ineqlin(arange(n1 + 1,dot(2,n1)))
# .\Num_Stab_Approx.m:218
                                        # are given by the difference phi_plus - phi_minus
    
    # 10. Infer the regression coefficients in the original regression with
# unnormalized data
#----------------------------------------------------------------------
    if (normalize == 1) or (RM >= 5):
        B[arange(2,n),arange()]=multiply(dot((1.0 / std(X(arange(),arange(2,n))).T),std(Y)),B)
# .\Num_Stab_Approx.m:230
        # of the intercept
        B[1,arange()]=mean(Y) - dot(mean(X(arange(),arange(2,n))),B(arange(2,n),arange()))
# .\Num_Stab_Approx.m:233
    