# Generated with SMOP  0.41
from libsmop import *
# .\deaadditprofit.m

    
@function
def deaadditprofit(X=None,Y=None,varargin=None,*args,**kwargs):
    varargin = deaadditprofit.varargin
    nargin = deaadditprofit.nargin

    #DEAALLOC Data envelopment analysis additive profit inefficiency
#   Computes data envelopment analysis additive profit inefficiency
    
    #   out = DEAADDITPROFIT(X, Y, Name, Value) computes data envelopment analysis 
#   additive profit inefficiency model with inputs X and outputs Y.
#   Model properties are specified using one or more Name ,Value pair 
#   arguments.
    
    #   Additional properties:
#   - 'Xprice': input prices.
#   - 'Yprice': output prices.
#   - 'rhoX': input slacks weights. Default is MIP: 1 ./ X
#   - 'rhoY': output slacks weights. Default is MIP: 1 ./ Y.
#   - 'names': DMU names.
    
    #   Example
#     
#      addprofit = deaadditprofit(X, Y, 'Xprice', W, 'Yprice', C);
    
    #   See also DEAOUT, DEA, DEAADDIT, DEAALLOC
    
    #   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
#   http://www.deatoolbox.com
    
    #   Version: 1.0
#   LAST UPDATE: 26, April, 2017
    
    
    # Check size
    if size(X,1) != size(Y,1):
        error('Number of rows in X must be equal to number of rows in Y')
    
    
    # Get number of DMUs (n), inputs (m) and outputs (s)
    n,m=size(X,nargout=2)
# .\deaadditprofit.m:36
    s=size(Y,2)
# .\deaadditprofit.m:37
    
    options=getDEAoptions(n,varargin[arange()])
# .\deaadditprofit.m:40
    
    rts=options.rts
# .\deaadditprofit.m:43
    if 'crs' == (rts):
        error('Only \'vrs\' are allowed')
    
    
    # If evaluate DMU at different X or Y
    if logical_not(isempty(options.Xeval)):
        Xeval=options.Xeval
# .\deaadditprofit.m:53
    else:
        Xeval=copy(X)
# .\deaadditprofit.m:55
    
    
    if logical_not(isempty(options.Yeval)):
        Yeval=options.Yeval
# .\deaadditprofit.m:59
    else:
        Yeval=copy(Y)
# .\deaadditprofit.m:61
    
    
    if size(Xeval,1) != size(Yeval,1):
        # Check size: rows
        error('Number of rows in Xref and Yref must be equal')
    
    
    if size(Xeval,2) != size(X,2):
        # Check columns Xref
        error('Number of columns in Xref and X must be equal')
    
    
    if size(Yeval,2) != size(Y,2):
        # Check columns Yref
        error('Number of columns in Yref and Y must be equal')
    
    
    neval=size(Xeval,1)
# .\deaadditprofit.m:79
    
    W=options.Xprice
# .\deaadditprofit.m:82
    P=options.Yprice
# .\deaadditprofit.m:83
    
    if logical_not(isempty(W)) and isempty(P):
        error('Both \'Xprice\' and \'Yprice\' must be specified')
    
    
    if isempty(W) and logical_not(isempty(P)):
        error('Both \'Xprice\' and \'Yprice\' must be specified')
    
    
    if logical_not(isempty(W)) and logical_not(isempty(P)):
        dispstr='names/X/Y/eff.T/eff.A/eff.P'
# .\deaadditprofit.m:95
    
    
    # Expand W and P if needed (if all firms have same prices and costs)
    if logical_not(isempty(W)) and size(W,1) == 1:
        W=repelem(W,neval,1)
# .\deaadditprofit.m:100
    
    
    if logical_not(isempty(P)) and size(P,1) == 1:
        P=repelem(P,neval,1)
# .\deaadditprofit.m:104
    
    
    # OPTIMIZATION OPTIONS:
    optimopts=options.optimopts
# .\deaadditprofit.m:108
    
    rhoX=options.rhoX
# .\deaadditprofit.m:111
    rhoY=options.rhoY
# .\deaadditprofit.m:112
    if isempty(rhoX):
        # rhoX = ones(size(X));
        # MIP
        rhoX=1 / X
# .\deaadditprofit.m:117
    
    if isempty(rhoY):
        # rhoY = ones(size(Y));
        # MIP
        rhoY=1 / Y
# .\deaadditprofit.m:122
    
    
    # TECHNICAL efficiency
    tech=deaaddit(X,Y,varargin[arange()],'rts','vrs','rhoX',rhoX,'rhoY',rhoY)
# .\deaadditprofit.m:127
    eff.T = copy(tech.eff)
# .\deaadditprofit.m:128
    
    
    # Create variables to store results
    lambda_=nan(neval,n)
# .\deaadditprofit.m:134
    Xeff=nan(neval,m)
# .\deaadditprofit.m:135
    Yeff=nan(neval,s)
# .\deaadditprofit.m:136
    
    for j in arange(1,neval).reshape(-1):
        # Objective function
        f=- concat([zeros(1,n),- W(j,arange()),P(j,arange())])
# .\deaadditprofit.m:142
        A=concat([[X.T,- eye(m,m),zeros(m,s)],[- Y.T,zeros(s,m),eye(s,s)]])
# .\deaadditprofit.m:145
        b=concat([[zeros(m,1)],[- zeros(s,1)]])
# .\deaadditprofit.m:147
        Aeq=concat([ones(1,n),zeros(1,m),zeros(1,s)])
# .\deaadditprofit.m:149
        beq=1
# .\deaadditprofit.m:150
        lb=zeros(1,n + m + s)
# .\deaadditprofit.m:151
        z=linprog(f,A,b,Aeq,beq,lb,[],[],optimopts)
# .\deaadditprofit.m:154
        lambda_[j,arange()]=z(arange(1,n))
# .\deaadditprofit.m:157
        Xeff[j,arange()]=z(arange(n + 1,n + m))
# .\deaadditprofit.m:158
        Yeff[j,arange()]=z(arange(n + m + 1,end()))
# .\deaadditprofit.m:159
    
    
    numerator=((sum(multiply(Yeff,P),2) - sum(multiply(Xeff,W),2)) - (sum(multiply(Y,P),2) - sum(multiply(X,W),2)))
# .\deaadditprofit.m:163
    denominator=min(min(concat([W / rhoX,P / rhoY])))
# .\deaadditprofit.m:164
    eff.P = copy(numerator / denominator)
# .\deaadditprofit.m:165
    
    eff.A = copy(eff.P - eff.T)
# .\deaadditprofit.m:168
    
    slack.X = copy(NaN)
# .\deaadditprofit.m:171
    slack.Y = copy(NaN)
# .\deaadditprofit.m:172
    Eflag=copy(NaN)
# .\deaadditprofit.m:174
    
    out=deaout('n',n,'neval',neval.T,'s',s,'m',m,'X',X,'Y',Y,'names',options.names,'model','additive-profit','orient','none','rts',rts,'lambda',lambda_,'slack',slack,'eff',eff,'Xeff',Xeff,'Yeff',Yeff,'exitflag',Eflag,'dispstr',dispstr,'Xprice',W,'Yprice',P)
# .\deaadditprofit.m:177
    return out
    
if __name__ == '__main__':
    pass
    