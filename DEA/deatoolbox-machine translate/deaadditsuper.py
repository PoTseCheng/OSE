# Generated with SMOP  0.41
from libsmop import *
# .\deaadditsuper.m

    
@function
def deaadditsuper(X=None,Y=None,varargin=None,*args,**kwargs):
    varargin = deaadditsuper.varargin
    nargin = deaadditsuper.nargin

    #DEAADDITSUPER Data envelopment analysis super efficiency additive model
#   Computes data envelopment analysis super efficiency additive model
    
    #   out = DEAADDITSUPER(X, Y, Name, Value) computes data envelopment analysis 
#   super efficiency additive model with inputs X and outputs Y. Model 
#   properties are specified using one or more Name ,Value pair arguments.
    
    #   Additional properties:
#   - 'rts': returns to sacle. Constant returns to scale 'crs', variable
#   returns to sacle 'vrs'.
#   - 'rhoX': input slacks weights. Default is MIP: 1 ./ X.
#   - 'rhoY': output slacks weights. Default is MIP: 1 ./ Y.
#   - 'names': DMU names.
    
    #   Example
#     
#      additsuper = deaadditsuper(X, Y, 'rts', 'vrs');
    
    #   See also DEAOUT, DEA, DEASCALE, DEAMALM, DEAADDIT, DEASUPER
    
    #   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
#   http://www.deatoolbox.com
    
    #   Version: 1.0
#   LAST UPDATE: 26, April, 2017
    
    # Check size
    if size(X,1) != size(Y,1):
        error('Number of rows in X must be equal to number of rows in Y')
    
    
    # TODO: Additive model (other function)
    
    # Get number of DMUs (n), inputs (m) and outputs (s)
    n,m=size(X,nargout=2)
# .\deaadditsuper.m:37
    s=size(Y,2)
# .\deaadditsuper.m:38
    neval=copy(n)
# .\deaadditsuper.m:39
    
    options=getDEAoptions(n,varargin[arange()])
# .\deaadditsuper.m:42
    orient=options.orient
# .\deaadditsuper.m:43
    
    if logical_not(strcmp(orient,'none')):
        error('Additive super-efficiency model orientation msut be none.')
    
    
    # Xeval, X and Yeval, Y must be equal in this function
    if logical_not(isempty(options.Xeval)) and size(options.Xeval) != size(X):
        error('Xeval and X must be equal')
    
    
    if logical_not(isempty(options.Yeval)) and size(options.Yeval) != size(Y):
        error('Yeval and Y must be equal')
    
    
    # OPTIMIZATION OPTIONS:
    optimopts=options.optimopts
# .\deaadditsuper.m:60
    
    lambda_=nan(neval,n - 1)
# .\deaadditsuper.m:63
    slackX=nan(neval,m)
# .\deaadditsuper.m:64
    slackY=nan(neval,s)
# .\deaadditsuper.m:65
    supereff=nan(n,1)
# .\deaadditsuper.m:66
    Xeff=nan(neval,m)
# .\deaadditsuper.m:67
    Yeff=nan(neval,s)
# .\deaadditsuper.m:68
    Eflag=nan(neval,1)
# .\deaadditsuper.m:69
    
    rts=options.rts
# .\deaadditsuper.m:72
    if 'crs' == (rts):
        #AeqRTS1 = [];
            #beqRTS1 = [];
        AeqRTS2super=[]
# .\deaadditsuper.m:78
        beqRTS2super=[]
# .\deaadditsuper.m:79
    else:
        if 'vrs' == (rts):
            #AeqRTS1 = [ones(1,n), 0];
            #beqRTS1 = 1;
            AeqRTS2super=concat([ones(1,n - 1),zeros(1,m),zeros(1,s)])
# .\deaadditsuper.m:84
            beqRTS2super=1
# .\deaadditsuper.m:85
    
    
    # SLACKS WEIGHTS
    rhoX=options.rhoX
# .\deaadditsuper.m:89
    rhoY=options.rhoY
# .\deaadditsuper.m:90
    if isempty(rhoX):
        # rhoX = ones(size(X));
        # MIP
        rhoX=1 / X
# .\deaadditsuper.m:95
    
    if isempty(rhoY):
        # rhoY = ones(size(Y));
        # MIP
        rhoY=1 / Y
# .\deaadditsuper.m:100
    
    
    # For each DMU
    for j in arange(1,n).reshape(-1):
        # ADDITIVE MODEL for each DMU
        tempdea=deaaddit(X,Y,varargin[arange()],'Xeval',X(j,arange()),'Yeval',Y(j,arange()))
# .\deaadditsuper.m:107
        if tempdea.eff < 1e-05:
            # ADDITIVE SUPER-EFFICIENCY
            # Objective Function
            fsuper=concat([zeros(1,n - 1),multiply(rhoX(j,arange()),ones(1,m)),multiply(rhoY(j,arange()),ones(1,s))])
# .\deaadditsuper.m:117
            lbsuper=zeros(n + m + s - 1,1)
# .\deaadditsuper.m:120
            others=arange(1,n)
# .\deaadditsuper.m:123
            others=others(others != j)
# .\deaadditsuper.m:124
            Asuper=concat([[X(others,arange()).T,- eye(m,m),zeros(m,s)],[- Y(others,arange()).T,zeros(s,m),- eye(s,s)]])
# .\deaadditsuper.m:127
            bsuper=concat([[X(j,arange()).T],[- Y(j,arange()).T]])
# .\deaadditsuper.m:129
            zsuper,__,exitflag=linprog(fsuper,Asuper,bsuper,AeqRTS2super,beqRTS2super,lbsuper,[],[],optimopts,nargout=3)
# .\deaadditsuper.m:130
            if exitflag != 1:
                if options.warning:
                    warning('Optimization exit flag: %i',exitflag)
            if isempty(zsuper):
                if options.warning:
                    warning('Optimization doesn\'t return a result. Results set to NaN.')
                zsuper=nan(n + m + s - 1,1)
# .\deaadditsuper.m:140
            lambda_[j,arange()]=zsuper(arange(1,n - 1))
# .\deaadditsuper.m:143
            slackX[j,arange()]=zsuper(arange(n,n + m - 1))
# .\deaadditsuper.m:144
            slackY[j,arange()]=zsuper(arange(n + m,n + m + s - 1))
# .\deaadditsuper.m:145
            supereff[j,arange()]=sum(multiply(rhoX(j,arange()),slackX(j,arange()))) + sum(multiply(rhoY(j,arange()),slackY(j,arange())))
# .\deaadditsuper.m:146
            Xeff[j,arange()]=NaN
# .\deaadditsuper.m:148
            Yeff[j,arange()]=NaN
# .\deaadditsuper.m:149
            Eflag[j]=exitflag
# .\deaadditsuper.m:151
        else:
            supereff[j]=NaN
# .\deaadditsuper.m:154
            lambda_[j,arange()]=NaN
# .\deaadditsuper.m:155
            slackX[j,arange()]=NaN
# .\deaadditsuper.m:156
            slackY[j,arange()]=NaN
# .\deaadditsuper.m:157
            Xeff[j,arange()]=NaN
# .\deaadditsuper.m:158
            Yeff[j,arange()]=NaN
# .\deaadditsuper.m:159
            Eflag[j]=tempdea.exitflag
# .\deaadditsuper.m:160
    
    
    # Slacks structure
    slack.X = copy(slackX)
# .\deaadditsuper.m:168
    slack.Y = copy(slackY)
# .\deaadditsuper.m:169
    
    out=deaout('n',n,'neval',neval.T,'s',s,'m',m,'X',X,'Y',Y,'names',options.names,'model','additive-supereff','orient',orient,'rts',rts,'lambda',lambda_,'slack',slack,'eff',supereff,'Xeff',Xeff,'Yeff',Yeff,'exitflag',Eflag,'dispstr','names/X/Y/slack.X/slack.Y/eff')
# .\deaadditsuper.m:173
    return out
    
if __name__ == '__main__':
    pass
    