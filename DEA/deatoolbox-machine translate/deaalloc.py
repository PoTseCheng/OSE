# Generated with SMOP  0.41
from libsmop import *
# .\deaalloc.m

    
@function
def deaalloc(X=None,Y=None,varargin=None,*args,**kwargs):
    varargin = deaalloc.varargin
    nargin = deaalloc.nargin

    #DEAALLOC Data envelopment analysis allocative model
#   Computes data envelopment analysis allocative model: cost, revenue and
#   profit.
    
    #   out = DEAALLOC(X, Y, Name, Value) computes data envelopment analysis 
#   allocative model (cost, revenue and profit) with inputs X and outputs Y.
#   Model properties are specified using one or more Name ,Value pair 
#   arguments.
    
    #   Additional properties:
#   - 'Xprice': input prices.
#   - 'Yprice': output prices.
#   - 'rts': returns to scale. Constant returns to scale 'crs', variable
#   returns to scale 'vrs'.
#   - 'Gx': input directions for profit model. Default is X.
#   - 'Gy': output directions for profit model. Default is Y.
#   - 'names': DMU names.
    
    #   Example
#     
#      cost = deaalloc(X, Y, 'Xprice', W);
#      revenuew = deaalloc(X, Y, 'Yprice', P);
#      profit = deaalloc(X, Y, 'Xprice', W, 'Yprice', C);
    
    #   See also DEAOUT, DEA, DEASCALE, DEAMALM, DEAADDIT, DEASUPER
    
    #   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
#   http://www.deatoolbox.com
    
    #   Version: 1.0
#   LAST UPDATE: 6, May, 2018
    
    
    # Check size
    if size(X,1) != size(Y,1):
        error('Number of rows in X must be equal to number of rows in Y')
    
    
    # Get number of DMUs (n), inputs (m) and outputs (s)
    n,m=size(X,nargout=2)
# .\deaalloc.m:41
    s=size(Y,2)
# .\deaalloc.m:42
    
    options=getDEAoptions(n,varargin[arange()])
# .\deaalloc.m:45
    
    rts=options.rts
# .\deaalloc.m:48
    
    if logical_not(isempty(options.Xeval)):
        Xeval=options.Xeval
# .\deaalloc.m:52
    else:
        Xeval=copy(X)
# .\deaalloc.m:54
    
    
    if logical_not(isempty(options.Yeval)):
        Yeval=options.Yeval
# .\deaalloc.m:58
    else:
        Yeval=copy(Y)
# .\deaalloc.m:60
    
    
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
# .\deaalloc.m:78
    
    W=options.Xprice
# .\deaalloc.m:81
    P=options.Yprice
# .\deaalloc.m:82
    
    if logical_not(isempty(W)) and isempty(P):
        model='allocative-cost'
# .\deaalloc.m:86
        orient='io'
# .\deaalloc.m:87
        dispstr='names/X/Xprice/Y/eff.T/eff.A/eff.C'
# .\deaalloc.m:88
    
    
    if isempty(W) and logical_not(isempty(P)):
        model='allocative-revenue'
# .\deaalloc.m:92
        orient='oo'
# .\deaalloc.m:93
        dispstr='names/X/Y/Yprice/eff.T/eff.A/eff.R'
# .\deaalloc.m:94
    
    
    if logical_not(isempty(W)) and logical_not(isempty(P)):
        model='allocative-profit'
# .\deaalloc.m:98
        orient='ddf'
# .\deaalloc.m:99
        dispstr='names/X/Xprice/Y/Yprice/eff.T/eff.A/eff.P'
# .\deaalloc.m:100
    
    
    # Error if RTS is not VRS for profit model
    if strcmp(model,'allocative-profit') and logical_not(strcmp(rts,'vrs')):
        error('RTS must be VRS for allocative profit model')
    
    
    # Expand W and P if needed (if all firms have same prices and costs)
    if logical_not(isempty(W)) and size(W,1) == 1:
        W=repelem(W,neval,1)
# .\deaalloc.m:110
    
    
    if logical_not(isempty(P)) and size(P,1) == 1:
        P=repelem(P,neval,1)
# .\deaalloc.m:114
    
    
    # OPTIMIZATION OPTIONS:
    optimopts=options.optimopts
# .\deaalloc.m:118
    
    lambda_=nan(neval,n)
# .\deaalloc.m:121
    Xeff=nan(neval,m)
# .\deaalloc.m:122
    Yeff=nan(neval,s)
# .\deaalloc.m:123
    
    if 'allocative-cost' == model:
        if 'crs' == (rts):
            AeqRTS1=[]
# .\deaalloc.m:133
            beqRTS1=[]
# .\deaalloc.m:134
        else:
            if 'vrs' == (rts):
                AeqRTS1=concat([ones(1,n),zeros(1,m)])
# .\deaalloc.m:136
                beqRTS1=1
# .\deaalloc.m:137
        # For each DMU
        for j in arange(1,neval).reshape(-1):
            # Objective function
            f=concat([zeros(1,n),W(j,arange())])
# .\deaalloc.m:144
            A=concat([[X.T,- eye(m,m)],[- Y.T,zeros(s,m)]])
# .\deaalloc.m:147
            b=concat([[zeros(m,1)],[- Yeval(j,arange()).T]])
# .\deaalloc.m:149
            Aeq=copy(AeqRTS1)
# .\deaalloc.m:151
            beq=copy(beqRTS1)
# .\deaalloc.m:152
            lb=zeros(1,n + m)
# .\deaalloc.m:153
            z=linprog(f,A,b,Aeq,beq,lb,[],[],optimopts)
# .\deaalloc.m:156
            lambda_[j,arange()]=z(arange(1,n))
# .\deaalloc.m:159
            Xeff[j,arange()]=z(arange(n + 1,end()))
# .\deaalloc.m:160
        # Cost efficiency
        eff.C = copy(sum(multiply(Xeff,W),2) / sum(multiply(X,W),2))
# .\deaalloc.m:165
        tempdea=dea(X,Y,varargin[arange()],'orient',orient)
# .\deaalloc.m:168
        eff.T = copy(tempdea.eff)
# .\deaalloc.m:169
        eff.A = copy(eff.C / eff.T)
# .\deaalloc.m:172
    else:
        if 'allocative-revenue' == model:
            if 'crs' == (rts):
                AeqRTS1=[]
# .\deaalloc.m:179
                beqRTS1=[]
# .\deaalloc.m:180
            else:
                if 'vrs' == (rts):
                    AeqRTS1=concat([ones(1,n),zeros(1,s)])
# .\deaalloc.m:182
                    beqRTS1=1
# .\deaalloc.m:183
            # For each DMU
            for j in arange(1,neval).reshape(-1):
                # Objective function
                f=- concat([zeros(1,n),P(j,arange())])
# .\deaalloc.m:190
                A=concat([X.T,- zeros(m,s),- Y.T,eye(s,s)])
# .\deaalloc.m:193
                b=concat([[Xeval(j,arange()).T],[- zeros(s,1)]])
# .\deaalloc.m:195
                Aeq=copy(AeqRTS1)
# .\deaalloc.m:197
                beq=copy(beqRTS1)
# .\deaalloc.m:198
                lb=zeros(1,n + s)
# .\deaalloc.m:199
                z=linprog(f,A,b,Aeq,beq,lb,[],[],optimopts)
# .\deaalloc.m:202
                lambda_[j,arange()]=z(arange(1,n))
# .\deaalloc.m:205
                Yeff[j,arange()]=z(arange(n + 1,end()))
# .\deaalloc.m:206
            # Revenue efficiency
            eff.R = copy(sum(multiply(Y,P),2) / sum(multiply(Yeff,P),2))
# .\deaalloc.m:211
            tempdea=dea(X,Y,varargin[arange()],'orient',orient)
# .\deaalloc.m:214
            eff.T = copy(1 / tempdea.eff)
# .\deaalloc.m:215
            eff.A = copy(eff.R / eff.T)
# .\deaalloc.m:218
        else:
            if 'allocative-profit' == model:
                # For each DMU
                for j in arange(1,neval).reshape(-1):
                    # Objective function
                    f=- concat([zeros(1,n),- W(j,arange()),P(j,arange())])
# .\deaalloc.m:227
                    A=concat([[X.T,- eye(m,m),zeros(m,s)],[- Y.T,zeros(s,m),eye(s,s)]])
# .\deaalloc.m:230
                    #     -Yeval(j,:)'];
                    b=concat([[zeros(m,1)],[- zeros(s,1)]])
# .\deaalloc.m:234
                    #beq = beqRTS1;
                    Aeq=concat([ones(1,n),zeros(1,m),zeros(1,s)])
# .\deaalloc.m:238
                    beq=1
# .\deaalloc.m:239
                    lb=zeros(1,n + m + s)
# .\deaalloc.m:240
                    z=linprog(f,A,b,Aeq,beq,lb,[],[],optimopts)
# .\deaalloc.m:243
                    lambda_[j,arange()]=z(arange(1,n))
# .\deaalloc.m:246
                    Xeff[j,arange()]=z(arange(n + 1,n + m))
# .\deaalloc.m:247
                    Yeff[j,arange()]=z(arange(n + m + 1,end()))
# .\deaalloc.m:248
                # Get directions
                Gx=options.Gx
# .\deaalloc.m:253
                Gy=options.Gy
# .\deaalloc.m:254
                if length(Gx) == 1:
                    Gx=repmat(Gx,size(X,1),size(X,2))
# .\deaalloc.m:257
                if length(Gy) == 1:
                    Gy=repmat(Gy,size(Y,1),size(Y,2))
# .\deaalloc.m:261
                if isempty(Gx):
                    Gx=copy(X)
# .\deaalloc.m:265
                if isempty(Gy):
                    Gy=copy(Y)
# .\deaalloc.m:269
                # Profit efficiency
            # Cambiado signo + por -
                eff.P = copy(((sum(multiply(Yeff,P),2) - sum(multiply(Xeff,W),2)) - (sum(multiply(Y,P),2) - sum(multiply(X,W),2))) / (sum(multiply(P,Gy),2) + sum(multiply(W,Gx),2)))
# .\deaalloc.m:274
                # Technical efficiency. DDF DEA model under VRS.
                tempdea=dea(X,Y,varargin[arange()],'orient',orient,'rts','vrs','Gx',Gx,'Gy',Gy)
# .\deaalloc.m:278
                eff.T = copy(tempdea.eff)
# .\deaalloc.m:279
                eff.A = copy(eff.P - eff.T)
# .\deaalloc.m:282
                # TODO: lambdas TEech or Profit???
    
    
    # Slacks structure
    slack.X = copy(NaN)
# .\deaalloc.m:290
    slack.Y = copy(NaN)
# .\deaalloc.m:291
    Eflag=copy(NaN)
# .\deaalloc.m:293
    
    out=deaout('n',n,'neval',neval.T,'s',s,'m',m,'X',X,'Y',Y,'names',options.names,'model',model,'orient',orient,'rts',rts,'lambda',lambda_,'slack',slack,'eff',eff,'Xeff',Xeff,'Yeff',Yeff,'exitflag',Eflag,'dispstr',dispstr,'Xprice',W,'Yprice',P)
# .\deaalloc.m:296
    return out
    
if __name__ == '__main__':
    pass
    