# Generated with SMOP  0.41
from libsmop import *
# .\deasuper.m

    
@function
def deasuper(X=None,Y=None,varargin=None,*args,**kwargs):
    varargin = deasuper.varargin
    nargin = deasuper.nargin

    #DEASUPER Data envelopment analysis super efficiency radial and directional
#   Computes data envelopment analysis super efficiency radial and 
#   directional model
    
    #   out = DEASUPER(X, Y, Name, Value) computes data envelopment analysis 
#   super efficiency model with inputs X and outputs Y. Model properties 
#   are specified using one or more Name ,Value pair arguments.
    
    #   Additional properties:
#   - 'orient': orientation. Input oriented 'io', output oriented 'oo', 
#   directional distane function 'ddf'.
#   - 'rts': returns to sacle. Constant returns to scale 'crs', variable
#   returns to sacle 'vrs'.
#   - 'Gx': input directions for 'ddf' orientation. Default is X.
#   - 'Gy': output directions for 'ddf' orientation. Default is Y.
#   - 'names': DMU names.
    
    #   Example
#     
#      iosuper = deasuper(X, Y, 'orient', 'io');
    
    #   See also DEAOUT, DEA, DEASCALE, DEAMALM, DEAADDIT, DEAADDITSUPER
    
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
# .\deasuper.m:40
    s=size(Y,2)
# .\deasuper.m:41
    neval=copy(n)
# .\deasuper.m:42
    
    options=getDEAoptions(n,varargin[arange()])
# .\deasuper.m:45
    orient=options.orient
# .\deasuper.m:46
    rts=options.rts
# .\deasuper.m:47
    
    if strcmp(orient,'none'):
        error('super-efficiency model must be oriented.')
    
    
    # Xeval, X and Yeval, Y must be equal in this function
    if logical_not(isempty(options.Xeval)) and size(options.Xeval) != size(X):
        error('Xeval and X must be equal')
    
    
    if logical_not(isempty(options.Yeval)) and size(options.Yeval) != size(Y):
        error('Yeval and Y must be equal')
    
    
    # If DDF
    if strcmp(orient,'ddf'):
        # Get directions
        Gx=options.Gx
# .\deasuper.m:66
        Gy=options.Gy
# .\deasuper.m:67
        if length(Gx) == 1:
            Gx=repmat(Gx,size(X,1),size(X,2))
# .\deasuper.m:70
        if length(Gy) == 1:
            Gy=repmat(Gy,size(Y,1),size(Y,2))
# .\deasuper.m:74
        if isempty(Gx):
            Gx=copy(X)
# .\deasuper.m:78
        if isempty(Gy):
            Gy=copy(Y)
# .\deasuper.m:82
    
    
    # Create variable to store results
    lambda_=nan(neval,n - 1)
# .\deasuper.m:88
    slackX=nan(neval,m)
# .\deasuper.m:89
    slackY=nan(neval,s)
# .\deasuper.m:90
    supereff=nan(n,1)
# .\deasuper.m:91
    Xeff=nan(neval,m)
# .\deasuper.m:92
    Yeff=nan(neval,s)
# .\deasuper.m:93
    Eflag=nan(neval,2)
# .\deasuper.m:94
    
    for j in arange(1,n).reshape(-1):
        # Evaluate each DMU w.r.t all without including itself
        others=arange(1,n)
# .\deasuper.m:100
        others=others(others != j)
# .\deasuper.m:101
        if strcmp(orient,'ddf'):
            # DDF
            tempdea=dea(X(others,arange()),Y(others,arange()),'orient',options.orient,'rts',options.rts,'Xeval',X(j,arange()),'Yeval',Y(j,arange()),'Gx',Gx(j,arange()),'Gy',Gy(j,arange()))
# .\deasuper.m:105
        else:
            tempdea=dea(X(others,arange()),Y(others,arange()),'orient',options.orient,'rts',options.rts,'Xeval',X(j,arange()),'Yeval',Y(j,arange()))
# .\deasuper.m:111
        supereff[j]=tempdea.eff(1)
# .\deasuper.m:119
        lambda_[j,arange()]=tempdea.lambda(1,arange())
# .\deasuper.m:120
        slackX[j,arange()]=tempdea.slack.X(1,arange())
# .\deasuper.m:121
        slackY[j,arange()]=tempdea.slack.Y(1,arange())
# .\deasuper.m:122
        Xeff[j,arange()]=tempdea.Xeff(1,arange())
# .\deasuper.m:123
        Yeff[j,arange()]=tempdea.Yeff(1,arange())
# .\deasuper.m:124
        Eflag[j,arange()]=tempdea.exitflag
# .\deasuper.m:125
    
    
    # Slacks structure
    slack.X = copy(slackX)
# .\deasuper.m:130
    slack.Y = copy(slackY)
# .\deasuper.m:131
    
    if strcmp(orient,'ddf'):
        model='directional-supereff'
# .\deasuper.m:135
    else:
        model='radial-supereff'
# .\deasuper.m:137
    
    
    # SAVE results and input data
    out=deaout('n',n,'neval',neval.T,'s',s,'m',m,'X',X,'Y',Y,'names',options.names,'model',model,'orient',orient,'rts',rts,'lambda',lambda_,'slack',slack,'eff',supereff,'Xeff',Xeff,'Yeff',Yeff,'exitflag',Eflag,'dispstr','names/X/Y/eff/slack.X/slack.Y')
# .\deasuper.m:141
    return out
    
if __name__ == '__main__':
    pass
    