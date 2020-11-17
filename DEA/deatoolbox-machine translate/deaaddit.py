# Generated with SMOP  0.41
from libsmop import *
# .\deaaddit.m

    
@function
def deaaddit(X=None,Y=None,varargin=None,*args,**kwargs):
    varargin = deaaddit.varargin
    nargin = deaaddit.nargin

    #DEAADDIT Data envelopment analysis weighted additive model
#   Computes data envelopment analysis weighted additive model
    
    #   out = DEAADDIT(X, Y, Name, Value) computes data envelopment analysis
#   weighted additive model with inputs X and outputs Y. Model properties 
#   are specified using one or more Name ,Value pair arguments. If weights
#   'rhoX' and 'rhoY' are not specified, the Measure of Inefficiency 
#   Proportions (MIP) program is computed.
    
    #   Additional properties:
#   - 'rts': returns to scale. Constant returns to scale 'crs', variable
#   returns to scale 'vrs'.
#   - 'rhoX': input slacks weights. Default is MIP: 1 ./ X.
#   - 'rhoY': output slacks weights. Default is MIP: 1 ./ Y.
#   - 'names': DMU names.
    
    #   Advanced parameters:
#   - 'Xeval: inputs to evaluate if different from X.
#   - 'Yeval': outputs to evaluate if different from Y.
    
    #   Example
#     
#      add = deaaddit(X, Y, 'rts', 'vrs');
    
    #   See also DEAOUT, DEA, DEAADDITSUPER
    
    #   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
#   http://www.deatoolbox.com
    
    #   Version: 1.0
#   LAST UPDATE: 26, April, 2017
    
    # Check size
    if size(X,1) != size(Y,1):
        error('Number of rows in X must be equal to number of rows in Y')
    
    
    # Get number of DMUs (n), inputs (m) and outputs (s)
    n,m=size(X,nargout=2)
# .\deaaddit.m:41
    s=size(Y,2)
# .\deaaddit.m:42
    
    options=getDEAoptions(n,varargin[arange()])
# .\deaaddit.m:45
    
    orient=options.orient
# .\deaaddit.m:48
    if logical_not(strcmp(options.orient,'none')):
        error('Additive model is non-oriented')
    
    
    # If evaluate DMU at different X or Y
    if logical_not(isempty(options.Xeval)):
        Xeval=options.Xeval
# .\deaaddit.m:55
    else:
        Xeval=copy(X)
# .\deaaddit.m:57
    
    
    if logical_not(isempty(options.Yeval)):
        Yeval=options.Yeval
# .\deaaddit.m:61
    else:
        Yeval=copy(Y)
# .\deaaddit.m:63
    
    
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
# .\deaaddit.m:81
    
    rts=options.rts
# .\deaaddit.m:84
    if 'crs' == (rts):
        #AeqRTS1 = [];
            #beqRTS1 = [];
        AeqRTS2=[]
# .\deaaddit.m:90
        beqRTS2=[]
# .\deaaddit.m:91
    else:
        if 'vrs' == (rts):
            #AeqRTS1 = [ones(1,n), 0];
            #beqRTS1 = 1;
            AeqRTS2=concat([ones(1,n),zeros(1,m),zeros(1,s)])
# .\deaaddit.m:96
            beqRTS2=1
# .\deaaddit.m:97
    
    
    # SLACKS WEIGHTS
    rhoX=options.rhoX
# .\deaaddit.m:101
    rhoY=options.rhoY
# .\deaaddit.m:102
    if isempty(rhoX):
        # rhoX = ones(size(X));
        # MIP
        rhoX=1 / X
# .\deaaddit.m:107
    
    if isempty(rhoY):
        # rhoY = ones(size(Y));
        # MIP
        rhoY=1 / Y
# .\deaaddit.m:112
    
    
    # OBJECTIVE FUNCTION
    # n zeros for \lambda
    # -ones for input slacks: m
    # -ones for output slacks: s    
    # Moved to the for loop to include weights
    
    # LOWER BOUNDS
    # Zero for n, m, s
    lb=zeros(n + m + s,1)
# .\deaaddit.m:123
    
    # Rows: n DMUs
    # Cols: n DMUs + m input slacks + s output slacks
    Z=zeros(neval,n + m + s)
# .\deaaddit.m:128
    Eflag=NaN(neval,1)
# .\deaaddit.m:129
    dualeqlin=nan(neval,m + s + logical_not(isempty(beqRTS2)))
# .\deaaddit.m:130
    
    optimopts=options.optimopts
# .\deaaddit.m:133
    
    # Solve linear problem for each DMU
    for j in arange(1,neval).reshape(-1):
        # Objective function with weights
        f=- concat([zeros(1,n),multiply(rhoX(j,arange()),ones(1,m)),multiply(rhoY(j,arange()),ones(1,s))])
# .\deaaddit.m:140
        Aeq=concat([[X.T,eye(m,m),zeros(m,s)],[Y.T,zeros(s,m),- eye(s,s)],[AeqRTS2]])
# .\deaaddit.m:142
        beq=concat([[Xeval(j,arange()).T],[Yeval(j,arange()).T],[beqRTS2]])
# .\deaaddit.m:145
        z,__,exitflag,__,dualz=linprog(f,[],[],Aeq,beq,lb,[],[],optimopts,nargout=5)
# .\deaaddit.m:146
        if exitflag != 1:
            if options.warning:
                warning('Optimization exit flag: %i',exitflag)
        if isempty(z):
            if options.warning:
                warning('Optimization doesn\'t return a result. Results set to NaN.')
            z=nan(n + m + s,1)
# .\deaaddit.m:156
            dualz.eqlin = copy(nan(1,m + s + logical_not(isempty(beqRTS2))))
# .\deaaddit.m:157
        Z[j,arange()]=z
# .\deaaddit.m:159
        Eflag[j]=exitflag
# .\deaaddit.m:160
        dualeqlin[j,arange()]=dualz.eqlin
# .\deaaddit.m:161
    
    
    # Get results
    lambda_=Z(arange(),arange(1,n))
# .\deaaddit.m:166
    slackX=Z(arange(),arange(n + 1,n + m))
# .\deaaddit.m:167
    slackY=Z(arange(),arange(n + m + 1,n + m + s))
# .\deaaddit.m:168
    eff=sum(multiply(rhoX(arange(1,neval),arange()),slackX),2) + sum(multiply(rhoY(arange(1,neval),arange()),slackY),2)
# .\deaaddit.m:169
    
    Xeff=Xeval - slackX
# .\deaaddit.m:172
    Yeff=Yeval + slackY
# .\deaaddit.m:173
    
    slack.X = copy(slackX)
# .\deaaddit.m:176
    slack.Y = copy(slackY)
# .\deaaddit.m:177
    
    dual.X = copy(dualeqlin(arange(),arange(1,m)))
# .\deaaddit.m:180
    dual.Y = copy(- dualeqlin(arange(),arange(m + 1,m + s)))
# .\deaaddit.m:181
    if logical_not(isempty(beqRTS2)):
        dual.rts = copy(- dualeqlin(arange(),arange(m + s + 1,m + s + 1)))
# .\deaaddit.m:183
    else:
        dual.rts = copy(nan(neval,1))
# .\deaaddit.m:185
    
    
    # SAVE results and input data
    out=deaout('n',n,'neval',neval.T,'s',s,'m',m,'X',X,'Y',Y,'names',options.names,'model','additive','orient',orient,'rts',rts,'lambda',lambda_,'slack',slack,'eff',eff,'Xeff',Xeff,'Yeff',Yeff,'dual',dual,'exitflag',Eflag,'dispstr','names/X/Y/slack.X/slack.Y/eff')
# .\deaaddit.m:189
    return out
    
if __name__ == '__main__':
    pass
    