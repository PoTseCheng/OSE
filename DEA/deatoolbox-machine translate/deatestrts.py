# Generated with SMOP  0.41
from libsmop import *
# .\deatestrts.m

    
@function
def deatestrts(X=None,Y=None,varargin=None,*args,**kwargs):
    varargin = deatestrts.varargin
    nargin = deatestrts.nargin

    #DEATESTRTS Data envelopment analysis test of returns to scale
#   Computes data envelopment analysis test of returns to scale (RTS)
    
    #   [S, SB, pvalue, critval] = DEATEST(X, Y, Name, Value) computes data 
#   envelopment analysis test of returns to scale with inputs X and outputs
#   Y. Model properties are specified using one or more Name ,Value pair 
#   arguments. The function returns the test statistic (S), the bootstrapped
#   statistic (SB), the p-value (pvalue) and the critical value (critval).
    
    #   Additional properties:
#   - 'orient': orientation. Input oriented 'io', output oriented 'oo'.
#   - 'names': DMU names.
#   - 'nreps': number of bootstrap replications. Default is 200.
#   - 'alpha': alpha value for confidence intervals. Default is 0.05.
#   - 'disp': set to 1 to display results on screen. Default is 0.
    
    #   Example
#     
#      [S, SB, pvalue, critval] = deatestrts(X, Y, 'orient', 'io');
    
    #   See also DEAOUT, DEA, DEASCALE, DEABOOT
    
    #   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
#   http://www.deatoolbox.com
    
    #   Version: 1.0
#   LAST UPDATE: 16, March, 2016
    
    # Get number of DMU's
    n=size(Y,1)
# .\deatestrts.m:32
    
    options=getDEAoptions(n,varargin[arange()])
# .\deatestrts.m:35
    nreps=options.nreps
# .\deatestrts.m:36
    alph=options.alpha
# .\deatestrts.m:37
    
    if logical_not(isempty(options.Xeval)) and size(options.Xeval) != size(X):
        error('Xeval and X must be equal')
    
    
    if logical_not(isempty(options.Yeval)) and size(options.Yeval) != size(Y):
        error('Yeval and Y must be equal')
    
    # Compute CRS DEA
    crs=dea(X,Y,varargin[arange()],'rts','crs')
# .\deatestrts.m:49
    crs_eff=crs.eff
# .\deatestrts.m:50
    
    vrs=dea(X,Y,varargin[arange()],'rts','vrs')
# .\deatestrts.m:53
    vrs_eff=vrs.eff
# .\deatestrts.m:54
    
    if cellarray(['io']) == (options.orient):
        S=sum(crs_eff) / sum(vrs_eff)
# .\deatestrts.m:59
    else:
        if cellarray(['oo']) == (options.orient):
            S=sum(1 / crs_eff) / sum(1 / vrs_eff)
# .\deatestrts.m:61
        else:
            if cellarray(['ddf']) == (options.orient):
                #scaleeff = crs_eff - vrs_eff;
                error('DEA rts test not available for \'ddf\'')
    
    # Bootstrap
    rset=rng()
# .\deatestrts.m:68
    
    crsB=deaboot(X,Y,varargin[arange()],'rts','crs','nreps',nreps,'alpha',alph)
# .\deatestrts.m:69
    crsB=crsB.eff.Boot
# .\deatestrts.m:70
    rng(rset)
    
    vrsB=deaboot(X,Y,varargin[arange()],'rts','vrs','nreps',nreps,'alpha',alph,'effRef',crs_eff)
# .\deatestrts.m:72
    vrsB=vrsB.eff.Boot
# .\deatestrts.m:73
    if strcmp(options.orient,'oo'):
        crsB=1 / crsB
# .\deatestrts.m:76
        vrsB=1 / vrsB
# .\deatestrts.m:77
    
    
    # Bootstrapped statistic
    SB=sum(crsB) / sum(vrsB)
# .\deatestrts.m:81
    
    lower=sum(SB < S)
# .\deatestrts.m:84
    pvalue=(lower + 1) / nreps
# .\deatestrts.m:85
    
    SBsorted=sort(SB)
# .\deatestrts.m:88
    critval=SBsorted(floor(dot(alph,nreps)))
# .\deatestrts.m:89
    
    if options.disp:
        fprintf('_______________________________\\n')
        fprintf('<strong>DEA Test of RTS</strong>\\n\\n')
        fprintf('H0: Globally CRS \\n')
        fprintf('H1: VRS \\n\\n')
        fprintf('Bootstrap replications: %i \\n',nreps)
        fprintf('Significance level: %4.2f \\n \\n',alph)
        fprintf('S statistic: %7.4f \\n',S)
        fprintf('Critical value: %7.4f \\n',critval)
        fprintf('p-value: %7.4f \\n',pvalue)
    
    
    return S,SB,pvalue,critval
    
if __name__ == '__main__':
    pass
    