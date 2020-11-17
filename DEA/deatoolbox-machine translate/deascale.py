# Generated with SMOP  0.41
from libsmop import *
# .\deascale.m

    
@function
def deascale(X=None,Y=None,varargin=None,*args,**kwargs):
    varargin = deascale.varargin
    nargin = deascale.nargin

    #DEASCALE Data envelopment analysis scale efficiency
#   Computes data envelopment analysis scale efficiency for radial and 
#   directional model
    
    #   out = DEASCALE(X, Y, Name, Value) computes data envelopment analysis 
#   scale efficiency model with inputs X and outputs Y. Model properties
#   are specified using one or more Name ,Value pair arguments.
    
    #   Additional properties:
#   - 'orient': orientation. Input oriented 'io', output oriented 'oo', 
#   directional distane function 'ddf'.
#   - 'Gx': input directions for 'ddf' orientation. Default is X.
#   - 'Gy': output directions for 'ddf' orientation. Default is Y.
#   - 'names': DMU names.
    
    #   Advanced parameters:
#   - 'Xeval: inputs to evaluate if different from X.
#   - 'Yeval': outputs to evaluate if different from Y.
    
    #   Example
#     
#      io_scale = deascale(X, Y, 'orient', 'io');
    
    #   See also DEAOUT, DEA, DEAMALM, DEAADDIT, DEASUPER
    
    #   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
#   http://www.deatoolbox.com
    
    #   Version: 1.0
#   LAST UPDATE: 9, May, 2017
    
    # Get number of DMU's
    n=size(Y,1)
# .\deascale.m:35
    
    options=getDEAoptions(n,varargin[arange()])
# .\deascale.m:38
    
    if logical_not(isempty(options.Xeval)):
        if size(options.Xeval) != size(X):
            error('Xeval and X must be of the same size')
    
    
    if logical_not(isempty(options.Yeval)):
        if size(options.Yeval) != size(Y):
            error('Yeval and Y must be of the same size')
    
    # Compute CRS DEA
    crs=dea(X,Y,varargin[arange()],'rts','crs')
# .\deascale.m:54
    crs_eff=crs.eff
# .\deascale.m:55
    
    vrs=dea(X,Y,varargin[arange()],'rts','vrs')
# .\deascale.m:58
    vrs_eff=vrs.eff
# .\deascale.m:59
    
    if cellarray(['io','oo']) == (options.orient):
        scaleeff=crs_eff / vrs_eff
# .\deascale.m:64
    else:
        if cellarray(['ddf']) == (options.orient):
            scaleeff=crs_eff - vrs_eff
# .\deascale.m:66
    
    
    # Efficiency
    eff.crs = copy(crs_eff)
# .\deascale.m:70
    eff.vrs = copy(vrs_eff)
# .\deascale.m:71
    eff.scale = copy(scaleeff)
# .\deascale.m:72
    
    neval=vrs.neval
# .\deascale.m:75
    s=vrs.s
# .\deascale.m:76
    m=vrs.m
# .\deascale.m:77
    slack.X = copy(NaN)
# .\deascale.m:78
    slack.Y = copy(NaN)
# .\deascale.m:79
    
    Eflag=nan(neval,4)
# .\deascale.m:82
    Eflag[arange(),arange(1,2)]=crs.exitflag
# .\deascale.m:83
    Eflag[arange(),arange(3,4)]=crs.exitflag
# .\deascale.m:84
    
    out=deaout('n',n,'neval',neval.T,'s',s,'m',m,'X',X,'Y',Y,'names',options.names,'model','radial','orient',options.orient,'rts','scaleeff','lambda',NaN,'slack',slack,'eff',eff,'Xeff',NaN,'Yeff',NaN,'exitflag',Eflag,'dispstr','names/eff.crs/eff.vrs/eff.scale')
# .\deascale.m:87
    return out
    
if __name__ == '__main__':
    pass
    