# Generated with SMOP  0.41
from libsmop import *
# .\dea.m

    
@function
def dea(X=None,Y=None,varargin=None,*args,**kwargs):
    varargin = dea.varargin
    nargin = dea.nargin

    #DEA Data envelopment analysis radial and directional model
#   Computes data envelopment analysis radial and directional model
    
    #   out = DEA(X, Y, Name, Value) computes data envelopment analysis model
#   with inputs X and outputs Y. Model properties are specified using 
#   one or more Name ,Value pair arguments.
    
    #   Additional properties:
#   - 'orient': orientation. Input oriented 'io', output oriented 'oo', 
#   directional distane function 'ddf'.
#   - 'rts': returns to scale. Constant returns to scale 'crs', variable
#   returns to scale 'vrs'.
#   - 'Gx': input directions for 'ddf' orientation. Default is Xeval.
#   - 'Gy': output directions for 'ddf' orientation. Default is Yeval.
#   - 'names': DMU names.
#   - 'secondstep': 1 to compute input and output slacks. Default is 1.
    
    #   Advanced parameters:
#   - 'Xeval: inputs to evaluate if different from X.
#   - 'Yeval': outputs to evaluate if different from Y.
    
    #   Example
#     
#      io = dea(X, Y, 'orient', 'io');
#      oo_vrs = dea(X, Y, 'orient', 'oo', 'rts', 'vrs');
#      ddf = dea(X, Y, 'ddf', 'Gx', X, 'Gy', Y);
    
    #   See also DEAOUT, DEASCALE, DEAMALM, DEAADDIT, DEASUPER
    
    #   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
#   http://www.deatoolbox.com
    
    #   Version: 1.0
#   LAST UPDATE: 1, August, 2017
    
    # Check size
    if size(X,1) != size(Y,1):
        error('Number of rows in X must be equal to number of rows in Y')
    
    
    # Get number of DMUs (n), inputs (m) and outputs (s)
    n,m=size(X,nargout=2)
# .\dea.m:44
    s=size(Y,2)
# .\dea.m:45
    
    options=getDEAoptions(n,varargin[arange()])
# .\dea.m:48
    
    if cellarray(['io','input']) == (options.orient):
        orient='io'
# .\dea.m:53
    else:
        if cellarray(['oo','output']) == (options.orient):
            orient='oo'
# .\dea.m:55
        else:
            if cellarray(['ddf']) == (options.orient):
                orient='ddf'
# .\dea.m:57
            else:
                if cellarray(['none']) == (options.orient):
                    error('Radial Model is Oriented')
    
    
    # RETURNS TO SCALE
    rts=options.rts
# .\dea.m:63
    if 'crs' == (rts):
        AeqRTS1=[]
# .\dea.m:66
        beqRTS1=[]
# .\dea.m:67
        AeqRTS2=[]
# .\dea.m:69
        beqRTS2=[]
# .\dea.m:70
    else:
        if 'vrs' == (rts):
            AeqRTS1=concat([ones(1,n),0])
# .\dea.m:72
            beqRTS1=1
# .\dea.m:73
            AeqRTS2=concat([ones(1,n),zeros(1,m),zeros(1,s)])
# .\dea.m:75
            beqRTS2=1
# .\dea.m:76
    
    
    # If evaluate DMU at different X or Y
    if logical_not(isempty(options.Xeval)):
        Xeval=options.Xeval
# .\dea.m:81
    else:
        Xeval=copy(X)
# .\dea.m:83
    
    
    if logical_not(isempty(options.Yeval)):
        Yeval=options.Yeval
# .\dea.m:87
    else:
        Yeval=copy(Y)
# .\dea.m:89
    
    
    if size(Xeval,1) != size(Yeval,1):
        # Check size: rows
        error('Number of rows in Xeval and Yeval must be equal')
    
    
    if size(Xeval,2) != size(X,2):
        # Check columns Xref
        error('Number of columns in Xeval and X must be equal')
    
    
    if size(Yeval,2) != size(Y,2):
        # Check columns Yref
        error('Number of columns in Yeval and Y must be equal')
    
    
    neval=size(Xeval,1)
# .\dea.m:107
    
    optimopts=options.optimopts
# .\dea.m:110
    
    lambda_=nan(neval,n)
# .\dea.m:113
    slackX=nan(neval,m)
# .\dea.m:114
    slackY=nan(neval,s)
# .\dea.m:115
    eff=nan(neval,1)
# .\dea.m:116
    Eflag=nan(neval,2)
# .\dea.m:117
    Xeff=nan(neval,m)
# .\dea.m:118
    Yeff=nan(neval,s)
# .\dea.m:119
    dualineqlin=nan(neval,m + s)
# .\dea.m:120
    dualeqlin=nan(neval,1)
# .\dea.m:121
    
    if 'io' == (orient):
        # For each DMU
        for j in arange(1,neval).reshape(-1):
            # FIRST STEP:
                # Objective function
            f=concat([zeros(1,n),1])
# .\dea.m:132
            A=concat([[X.T,- Xeval(j,arange()).T],[- Y.T,zeros(s,1)]])
# .\dea.m:135
            b=concat([[zeros(m,1)],[- Yeval(j,arange()).T]])
# .\dea.m:137
            Aeq=copy(AeqRTS1)
# .\dea.m:139
            beq=copy(beqRTS1)
# .\dea.m:140
            lb=zeros(1,n + 1)
# .\dea.m:141
            z,__,exitflag,__,dual=linprog(f,A,b,Aeq,beq,lb,[],[],optimopts,nargout=5)
# .\dea.m:144
            if exitflag != 1:
                if options.warning:
                    warning('DMU %i. First Step. Optimization exit flag: %i',j,exitflag)
            if isempty(z):
                if options.warning:
                    warning('DMU %i. First Step. Optimization doesn\'t return a result. Efficiency set to NaN.',j)
                z=nan(n + 1,1)
# .\dea.m:154
                dual.ineqlin = copy(nan(1,m + s))
# .\dea.m:155
                dual.eqlin = copy(nan(1,1))
# .\dea.m:156
            # Get efficiency
            theta=z(end())
# .\dea.m:160
            Eflag[j,1]=exitflag
# .\dea.m:161
            eff[j]=theta
# .\dea.m:162
            dualineqlin[j,arange()]=dual.ineqlin
# .\dea.m:165
            if logical_not(isempty(beqRTS1)):
                dualeqlin[j,1]=dual.eqlin
# .\dea.m:167
            else:
                dualeqlin[j,1]=NaN
# .\dea.m:169
            # SECOND STEP
            if (options.secondstep) and logical_not(isnan(theta)):
                # Objective function
                f=concat([zeros(1,n),- ones(1,m + s)])
# .\dea.m:177
                Aeq=concat([[X.T,eye(m,m),zeros(m,s)],[Y.T,zeros(s,m),- eye(s,s)],[AeqRTS2]])
# .\dea.m:180
                beq=concat([[multiply(theta,Xeval(j,arange()).T)],[Yeval(j,arange()).T],[beqRTS2]])
# .\dea.m:183
                lb=zeros(n + s + m,1)
# .\dea.m:186
                z,__,exitflag=linprog(f,[],[],Aeq,beq,lb,[],[],optimopts,nargout=3)
# .\dea.m:189
                if exitflag != 1:
                    if options.warning:
                        warning('DMU %i. Second Step. Optimization exit flag: %i',j,exitflag)
                if isempty(z):
                    if options.warning:
                        warning('DMU %i. Second Step. Optimization doesn\'t return a result. Results set to NaN.',j)
                    z=nan(n + m + s,1)
# .\dea.m:199
                # Get results
                lambda_[j,arange()]=z(arange(1,n))
# .\dea.m:203
                slackX[j,arange()]=z(arange(n + 1,n + m))
# .\dea.m:204
                slackY[j,arange()]=z(arange(n + m + 1,n + m + s))
# .\dea.m:205
                Eflag[j,2]=exitflag
# .\dea.m:206
                Xeff[j,arange()]=multiply(repmat(eff(j),1,m),Xeval(j,arange())) - slackX(j,arange())
# .\dea.m:209
                Yeff[j,arange()]=Yeval(j,arange()) + slackY(j,arange())
# .\dea.m:210
        # Compute efficient inputs and outputs
            # Xeff = repmat(eff, 1, m) .* Xeval - slackX;
            # Yeff = Yeval + slackY;
    else:
        if 'oo' == (orient):
            # For each DMU
            for j in arange(1,neval).reshape(-1):
                # FIRST STEP:
                # Objective function (maximize)
                f=- concat([zeros(1,n),1])
# .\dea.m:228
                A=concat([[X.T,zeros(m,1)],[- Y.T,Yeval(j,arange()).T]])
# .\dea.m:231
                b=concat([[Xeval(j,arange()).T],[zeros(s,1)]])
# .\dea.m:233
                Aeq=copy(AeqRTS1)
# .\dea.m:235
                beq=copy(beqRTS1)
# .\dea.m:236
                lb=zeros(1,n + 1)
# .\dea.m:237
                z,__,exitflag,__,dual=linprog(f,A,b,Aeq,beq,lb,[],[],optimopts,nargout=5)
# .\dea.m:240
                if exitflag != 1:
                    if options.warning:
                        warning('DMU %i. First Step. Optimization exit flag: %i',j,exitflag)
                if isempty(z):
                    if options.warning:
                        warning('DMU %i. First Step. Optimization doesn\'t return a result. Efficiency set to NaN.',j)
                    z=nan(n + 1,1)
# .\dea.m:250
                    dual.ineqlin = copy(nan(1,m + s))
# .\dea.m:251
                    dual.eqlin = copy(nan(1,1))
# .\dea.m:252
                # Get efficiency
                phi=z(end())
# .\dea.m:256
                eff[j]=phi
# .\dea.m:257
                Eflag[j,1]=exitflag
# .\dea.m:258
                dualineqlin[j,arange()]=dual.ineqlin
# .\dea.m:261
                if logical_not(isempty(beqRTS1)):
                    dualeqlin[j,1]=dual.eqlin
# .\dea.m:263
                else:
                    dualeqlin[j,1]=NaN
# .\dea.m:265
                # SECOND STEP
                if (options.secondstep) and logical_not(isnan(phi)):
                    # Objective function
                    f=- concat([zeros(1,n),ones(1,m + s)])
# .\dea.m:273
                    Aeq=concat([[X.T,eye(m,m),zeros(m,s)],[Y.T,zeros(s,m),- eye(s,s)],[AeqRTS2]])
# .\dea.m:276
                    beq=concat([[Xeval(j,arange()).T,multiply(phi,Yeval(j,arange()).T)],[beqRTS2]])
# .\dea.m:279
                    lb=zeros(n + s + m,1)
# .\dea.m:282
                    z,__,exitflag=linprog(f,[],[],Aeq,beq,lb,[],[],optimopts,nargout=3)
# .\dea.m:285
                    if exitflag != 1:
                        if options.warning:
                            warning('DMU %i. Second Step. Optimization exit flag: %i',j,exitflag)
                    if isempty(z):
                        if options.warning:
                            warning('DMU %i. Second Step. Optimization doesn\'t return a result. Results set to NaN.',j)
                        z=nan(n + m + s,1)
# .\dea.m:295
                    # Get results
                    lambda_[j,arange()]=z(arange(1,n))
# .\dea.m:299
                    slackX[j,arange()]=z(arange(n + 1,n + m))
# .\dea.m:300
                    slackY[j,arange()]=z(arange(n + m + 1,n + m + s))
# .\dea.m:301
                    Eflag[j,2]=exitflag
# .\dea.m:302
                    Xeff[j,arange()]=Xeval(j,arange()) - slackX(j,arange())
# .\dea.m:305
                    Yeff[j,arange()]=multiply(repmat(eff(j),1,s),Yeval(j,arange())) + slackY(j,arange())
# .\dea.m:306
            # Compute efficient inputs and outputs
            # Xeff = Xeval - slackX;
            # Yeff = repmat(eff, 1, s) .* Yeval + slackY;
        else:
            if 'ddf' == (orient):
                # Get directions
                Gx=options.Gx
# .\dea.m:320
                Gy=options.Gy
# .\dea.m:321
                if length(Gx) == 1:
                    Gx=repmat(Gx,size(X,1),size(X,2))
# .\dea.m:324
                else:
                    if size(Gx,1) == 1:
                        Gx=repmat(Gx,size(X,1),1)
# .\dea.m:326
                if length(Gy) == 1:
                    Gy=repmat(Gy,size(Y,1),size(Y,2))
# .\dea.m:330
                else:
                    if size(Gy,1) == 1:
                        Gy=repmat(Gy,size(Y,1),1)
# .\dea.m:332
                if isempty(Gx):
                    Gx=copy(Xeval)
# .\dea.m:336
                if isempty(Gy):
                    Gy=copy(Yeval)
# .\dea.m:340
                # For each DMU
                for j in arange(1,neval).reshape(-1):
                    # FIRST STEP:
                # Objective function (maximize)
                    f=- concat([zeros(1,n),1])
# .\dea.m:348
                    A=concat([[X.T,Gx(j,arange()).T],[- Y.T,Gy(j,arange()).T]])
# .\dea.m:351
                    b=concat([[Xeval(j,arange()).T],[- Yeval(j,arange()).T]])
# .\dea.m:353
                    Aeq=copy(AeqRTS1)
# .\dea.m:355
                    beq=copy(beqRTS1)
# .\dea.m:356
                    lb=concat([zeros(1,n),- inf])
# .\dea.m:358
                    z,__,exitflag,__,dual=linprog(f,A,b,Aeq,beq,lb,[],[],optimopts,nargout=5)
# .\dea.m:361
                    if exitflag != 1:
                        if options.warning:
                            warning('DMU %i. First Step. Optimization exit flag: %i',j,exitflag)
                    if isempty(z):
                        if options.warning:
                            warning('DMU %i. First Step. Optimization doesn\'t return a result. Efficiency set to NaN.',j)
                        z=nan(n + 1,1)
# .\dea.m:371
                        dual.ineqlin = copy(nan(1,m + s))
# .\dea.m:372
                        dual.eqlin = copy(nan(1,1))
# .\dea.m:373
                    # Get efficiency
                    beta=z(end())
# .\dea.m:377
                    eff[j]=beta
# .\dea.m:378
                    Eflag[j,1]=exitflag
# .\dea.m:379
                    dualineqlin[j,arange()]=dual.ineqlin
# .\dea.m:382
                    if logical_not(isempty(beqRTS1)):
                        dualeqlin[j,1]=dual.eqlin
# .\dea.m:384
                    else:
                        dualeqlin[j,1]=NaN
# .\dea.m:386
                    # SECOND STEP
                    if (options.secondstep) and logical_not(isnan(beta)):
                        # Objective function
                        f=- concat([zeros(1,n),ones(1,m + s)])
# .\dea.m:394
                        Aeq=concat([[X.T,eye(m,m),zeros(m,s)],[Y.T,zeros(s,m),- eye(s,s)],[AeqRTS2]])
# .\dea.m:397
                        beq=concat([[multiply(- beta,Gx(j,arange()).T) + Xeval(j,arange()).T,multiply(beta,Gy(j,arange()).T) + Yeval(j,arange()).T],[beqRTS2]])
# .\dea.m:400
                        lb=zeros(n + s + m,1)
# .\dea.m:403
                        z,__,exitflag=linprog(f,[],[],Aeq,beq,lb,[],[],optimopts,nargout=3)
# .\dea.m:406
                        if exitflag != 1:
                            if options.warning:
                                warning('DMU %i. Second Step. Optimization exit flag: %i',j,exitflag)
                        if isempty(z):
                            if options.warning:
                                warning('DMU %i. Second Step. Optimization doesn\'t return a result. Results set to NaN.',j)
                            z=nan(n + m + s,1)
# .\dea.m:416
                        # Get results
                        lambda_[j,arange()]=z(arange(1,n))
# .\dea.m:420
                        slackX[j,arange()]=z(arange(n + 1,n + m))
# .\dea.m:421
                        slackY[j,arange()]=z(arange(n + m + 1,n + m + s))
# .\dea.m:422
                        Eflag[j,2]=exitflag
# .\dea.m:423
                        Xeff[j,arange()]=Xeval(j,arange()) - multiply(repmat(eff(j),1,m),Gx(j,arange())) - slackX(j,arange())
# .\dea.m:426
                        Yeff[j,arange()]=Yeval(j,arange()) + multiply(repmat(eff(j),1,s),Gy(j,arange())) + slackY(j,arange())
# .\dea.m:427
    
    
    # Slacks structure
    slack.X = copy(slackX)
# .\dea.m:436
    slack.Y = copy(slackY)
# .\dea.m:437
    
    dual.X = copy(dualineqlin(arange(),arange(1,m)))
# .\dea.m:440
    dual.Y = copy(dualineqlin(arange(),arange(m + 1,m + s)))
# .\dea.m:441
    if logical_not(isempty(beqRTS2)):
        dual.rts = copy(ravel(dualeqlin))
# .\dea.m:443
    else:
        dual.rts = copy(nan(neval,1))
# .\dea.m:445
    
    
    # SAVE results and input data
    out=deaout('n',n,'neval',neval.T,'s',s,'m',m,'X',X,'Y',Y,'names',options.names,'model','radial','orient',orient,'rts',rts,'lambda',lambda_,'slack',slack,'eff',eff,'Xeff',Xeff,'Yeff',Yeff,'dual',dual,'exitflag',Eflag,'dispstr','names/X/Y/eff/slack.X/slack.Y')
# .\dea.m:449
    return out
    
if __name__ == '__main__':
    pass
    