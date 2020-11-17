# Generated with SMOP  0.41
from libsmop import *
# .\deaund.m

    
@function
def deaund(X=None,Y=None,Yu=None,varargin=None,*args,**kwargs):
    varargin = deaund.varargin
    nargin = deaund.nargin

    #DEAUND Data envelopment analysis with undesirable outputs.
#   Computes data envelopment analysis model with undesirable outputs.
    
    #   out = DEAUND(X, Y, Yu, Name, Value) computes data envelopment analysis 
#   model with inputs X, outputs Y, and undesirable outputs Yu. Model 
#   properties are specified using one or more Name ,Value pair arguments.
    
    #   Additional properties:
#   - 'names': DMU names.
#   - 'orient': orientation. Directional distane function with undesirable
#   outputs 'ddf' (Aparicio, Pastor and Zofio, 2013), default. Directional 
#   distance function with undesirable outputs 'ddf_cfg' (Chung, Fare and 
#   Grosskopf).
    
    #   Advanced parameters:
#   - 'Xeval: inputs to evaluate if different from X.
#   - 'Yeval': outputs to evaluate if different from Y.
    
    #   Example
#     
#      und = deaund(X, Y, Yu);
    
    #   See also DEAOUT, DEA, DEAMALMLUEN
    
    #   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
#   http://www.deatoolbox.com
    
    #   Version: 1.0
#   LAST UPDATE: 27, April, 2017
    
    # Check size
    if size(X,1) != size(Y,1):
        error('Number of rows in X must be equal to number of rows in Y')
    
    
    if size(Yu,1) != size(Y,1):
        error('Number of rows in Yu must be equal to number of rows in Y')
    
    
    # Get number of DMUs (n), inputs (m), outputs (s), and undesirable
    # outputs (r)
    n,m=size(X,nargout=2)
# .\deaund.m:44
    s=size(Y,2)
# .\deaund.m:45
    r=size(Yu,2)
# .\deaund.m:46
    
    
    # Get DEA options
    options=getDEAoptions(n,varargin[arange()])
# .\deaund.m:49
    
    if cellarray(['none']) == (options.orient):
        # Replace default 'none' orientation to 'ddf'
        orient='ddf'
# .\deaund.m:55
    else:
        if cellarray(['ddf']) == (options.orient):
            orient='ddf'
# .\deaund.m:57
        else:
            if cellarray(['ddf_cfg']) == (options.orient):
                orient='ddf_cfg'
# .\deaund.m:59
            else:
                error('Orientation for the undesarible outputs model must be ddf')
    
    
    # Distance functions
    if logical_not(isempty(options.Gx)) or logical_not(isempty(options.Gy)):
        error('Distance functions Gx, Gy, and Gyu are automatically assigned in undesirable outputs dea model.')
    
    
    # RETURNS TO SCALE
    rts=options.rts
# .\deaund.m:70
    if 'crs' == (rts):
        AeqRTS1=[]
# .\deaund.m:73
        beqRTS1=[]
# .\deaund.m:74
        AeqRTS2=[]
# .\deaund.m:76
        beqRTS2=[]
# .\deaund.m:77
    else:
        error('DEA model with undesirable outputs only available for crs')
    
    
    # If evaluate DMU at different X or Y
    if logical_not(isempty(options.Xeval)):
        Xeval=options.Xeval
# .\deaund.m:84
    else:
        Xeval=copy(X)
# .\deaund.m:86
    
    
    if logical_not(isempty(options.Yeval)):
        Yeval=options.Yeval
# .\deaund.m:90
    else:
        Yeval=copy(Y)
# .\deaund.m:92
    
    
    if logical_not(isempty(options.Yueval)):
        Yueval=options.Yueval
# .\deaund.m:96
    else:
        Yueval=copy(Yu)
# .\deaund.m:98
    
    
    if size(Xeval,1) != size(Yeval,1):
        # Check size: rows
        error('Number of rows in Xref and Yref must be equal')
    
    
    if size(Yueval,1) != size(Yeval,1):
        # Check size: rows
        error('Number of rows in Yuref and Yref must be equal')
    
    
    if size(Xeval,2) != size(X,2):
        # Check columns Xref
        error('Number of columns in Xref and X must be equal')
    
    
    if size(Yeval,2) != size(Y,2):
        # Check columns Yref
        error('Number of columns in Yref and Y must be equal')
    
    
    if size(Yueval,2) != size(Yu,2):
        # Check columns Yref
        error('Number of columns in Yuref and Yu must be equal')
    
    
    neval=size(Xeval,1)
# .\deaund.m:126
    
    optimopts=options.optimopts
# .\deaund.m:129
    
    lambda_=nan(neval,n)
# .\deaund.m:132
    slackX=nan(neval,m)
# .\deaund.m:133
    slackY=nan(neval,s)
# .\deaund.m:134
    slackYu=nan(neval,r)
# .\deaund.m:135
    eff=nan(neval,1)
# .\deaund.m:136
    Eflag=nan(neval,2)
# .\deaund.m:137
    Xeff=nan(neval,m)
# .\deaund.m:138
    Yeff=nan(neval,s)
# .\deaund.m:139
    Yueff=nan(neval,r)
# .\deaund.m:140
    
    if 'ddf' == (orient):
        # (Aparicio, Pastor and Zofio, 2013)
        # Get directions
            #G = options.ddfG;
            #H = options.ddfH;
        Gx=zeros(n,m)
# .\deaund.m:151
        Gy=copy(Yeval)
# .\deaund.m:152
        Gyu=copy(Yueval)
# .\deaund.m:153
        if length(Gx) == 1:
            Gx=repmat(Gx,size(X,1),size(X,2))
# .\deaund.m:156
        if length(Gy) == 1:
            Gy=repmat(Gy,size(Y,1),size(Y,2))
# .\deaund.m:160
        if length(Gyu) == 1:
            Gyu=repmat(Gyu,size(Y,1),size(Y,2))
# .\deaund.m:164
        maxYu=max(max(concat([[Yu],[Yueval]])))
# .\deaund.m:167
        for j in arange(1,neval).reshape(-1):
            # FIRST STEP:
                # Objective function (maximize)
            f=- concat([zeros(1,n),1])
# .\deaund.m:175
            A=concat([[X.T,Gx(j,arange()).T],[- Y.T,Gy(j,arange()).T],[Yu.T,Gyu(j,arange()).T],[zeros(r,n),- Yueval(j,arange()).T]])
# .\deaund.m:178
            b=concat([[Xeval(j,arange()).T],[- Yeval(j,arange()).T],[Yueval(j,arange()).T],[maxYu - Yueval(j,arange()).T]])
# .\deaund.m:182
            Aeq=copy(AeqRTS1)
# .\deaund.m:187
            beq=copy(beqRTS1)
# .\deaund.m:188
            lb=concat([zeros(1,n),- inf])
# .\deaund.m:190
            z,__,exitflag=linprog(f,A,b,Aeq,beq,lb,[],[],optimopts,nargout=3)
# .\deaund.m:194
            if exitflag != 1:
                if options.warning:
                    warning('DMU %i. First Step. Optimization exit flag: %i',j,exitflag)
            if isempty(z):
                if options.warning:
                    warning('DMU %i. First Step. Optimization doesn\'t return a result in First Step. Efficiency set to NaN.',j)
                z=nan(n + 1,1)
# .\deaund.m:204
            # Get efficiency
            beta=z(end())
# .\deaund.m:208
            eff[j]=beta
# .\deaund.m:209
            Eflag[j,1]=exitflag
# .\deaund.m:210
            if (options.secondstep) and logical_not(isnan(beta)):
                # Objective function
                f=- concat([zeros(1,n),ones(1,m + s + r)])
# .\deaund.m:216
                Aeq=concat([[X.T,eye(m,m),zeros(m,s),zeros(m,r)],[Y.T,zeros(s,m),- eye(s,s),zeros(s,r)],[Yu.T,zeros(r,m),zeros(r,s),eye(r,r)],[AeqRTS2]])
# .\deaund.m:219
                beq=concat([[multiply(- beta,Gx(j,arange()).T) + Xeval(j,arange()).T,multiply(beta,Gy(j,arange()).T) + Yeval(j,arange()).T],[multiply(- beta,Gyu(j,arange()).T) + Yueval(j,arange()).T],[beqRTS2]])
# .\deaund.m:223
                lb=zeros(n + s + m + r,1)
# .\deaund.m:227
                # Optimize
                z=linprog(f,[],[],Aeq,beq,lb,[],[],optimopts)
# .\deaund.m:231
                if exitflag != 1:
                    if options.warning:
                        warning('DMU %i. Second Step. Optimization exit flag: %i',j,exitflag)
                if isempty(z):
                    if options.warning:
                        warning('DMU %i. Second Step. Optimization doesn\'t return a result. Results set to NaN.',j)
                    z=nan(n + m + s + r,1)
# .\deaund.m:241
                # Get results
                lambda_[j,arange()]=z(arange(1,n))
# .\deaund.m:245
                slackX[j,arange()]=z(arange(n + 1,n + m))
# .\deaund.m:246
                slackY[j,arange()]=z(arange(n + m + 1,n + m + s))
# .\deaund.m:247
                slackYu[j,arange()]=z(arange(n + m + s + 1,n + m + s + r))
# .\deaund.m:248
                Eflag[j,2]=exitflag
# .\deaund.m:249
                Xeff[j,arange()]=Xeval(j,arange()) - multiply(repmat(eff(j),1,m),Gx(j,arange())) - slackX(j,arange())
# .\deaund.m:252
                Yeff[j,arange()]=Yeval(j,arange()) + multiply(repmat(eff(j),1,s),Gy(j,arange())) + slackY(j,arange())
# .\deaund.m:253
                Yueff[j,arange()]=Yueval(j,arange()) - multiply(repmat(eff(j),1,r),Gyu(j,arange())) - slackYu(j,arange())
# .\deaund.m:254
    else:
        if 'ddf_cfg' == (orient):
            # (Chung, Fare and Grosskopf)
            # For each DMU
            for j in arange(1,neval).reshape(-1):
                # FIRST STEP:
                # Objective function (maximize)
                f=- concat([zeros(1,n),1])
# .\deaund.m:267
                A=concat([[X.T,zeros(m,1)],[- Y.T,Yeval(j,arange()).T]])
# .\deaund.m:270
                b=concat([[Xeval(j,arange()).T],[- Yeval(j,arange()).T]])
# .\deaund.m:272
                Aeq=concat([Yu.T,Yueval(j,arange()).T])
# .\deaund.m:273
                beq=Yueval(j,arange()).T
# .\deaund.m:274
                lb=concat([zeros(1,n),- inf])
# .\deaund.m:275
                z,__,exitflag=linprog(f,A,b,Aeq,beq,lb,[],[],optimopts,nargout=3)
# .\deaund.m:279
                if exitflag != 1:
                    if options.warning:
                        warning('DMU %i. First Step. Optimization exit flag: %i',j,exitflag)
                if isempty(z):
                    if options.warning:
                        warning('DMU %i. First Step. Optimization doesn\'t return a result in First Step. Efficiency set to NaN.',j)
                    z=nan(n + 1,1)
# .\deaund.m:289
                # Get efficiency
                beta=z(end())
# .\deaund.m:293
                eff[j]=beta
# .\deaund.m:294
                Eflag[j,1]=exitflag
# .\deaund.m:295
                # Not available
                lambda_[j,arange()]=z(arange(1,n))
# .\deaund.m:299
                slackX[j,arange()]=nan(1,m)
# .\deaund.m:300
                slackY[j,arange()]=nan(1,s)
# .\deaund.m:301
                slackYu[j,arange()]=nan(1,r)
# .\deaund.m:302
                eff[j]=beta
# .\deaund.m:303
                Eflag[j,2]=nan(1,1)
# .\deaund.m:304
                Xeff[j,arange()]=nan(1,m)
# .\deaund.m:307
                Yeff[j,arange()]=nan(1,s)
# .\deaund.m:308
                Yueff[j,arange()]=nan(1,r)
# .\deaund.m:309
    
    
    # Slacks structure
    slack.X = copy(slackX)
# .\deaund.m:317
    slack.Y = copy(slackY)
# .\deaund.m:318
    slack.Yu = copy(slackYu)
# .\deaund.m:319
    
    out=deaout('n',n,'neval',neval.T,'s',s,'m',m,'X',X,'Y',Y,'names',options.names,'model','directional-undesirable','orient',orient,'rts',rts,'lambda',lambda_,'slack',slack,'eff',eff,'Xeff',Xeff,'Yeff',Yeff,'exitflag',Eflag,'dispstr','names/X/Y/Yu/eff/slack.X/slack.Y/slack.Yu','r',r,'Yu',Yu,'Yueff',Yueff)
# .\deaund.m:322
    return out
    
if __name__ == '__main__':
    pass
    