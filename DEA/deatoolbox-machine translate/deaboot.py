# Generated with SMOP  0.41
from libsmop import *
# .\deaboot.m

    
@function
def deaboot(X=None,Y=None,varargin=None,*args,**kwargs):
    varargin = deaboot.varargin
    nargin = deaboot.nargin

    #DEABOOT Data envelopment analysis bootstrap.
#   Computes data envelopment analysis bootstrap following Simar and Wilson
#   (1998) and Bogetoft and Otto (2001)
    
    #   out = DEABOOT(X, Y, Name, Value) computes data envelopment analysis
#   bootstrap model with inputs X and outputs Y. Model properties are 
#specified using one or more Name ,Value pair arguments.
    
    #   Additional properties:
#   - 'orient': orientation. Input oriented 'io', output oriented 'oo', 
#   directional distane function 'ddf'.
#   - 'rts': returns to sacle. Constant returns to scale 'crs', variable
#   returns to sacle 'vrs'.
#   - 'names': DMU names.
#   - 'nreps': number of bootstrap replications. Default is 200.
#   - 'alpha': alpha value for confidence intervals. Default is 0.05.
    
    #   Example
#     
#       io_b = deaboot(X, Y, 'orient', 'io', 'nreps', 200);
#       deadisp(io_b);
    
    #   See also DEAOUT, DEA
    
    #   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
#   http://www.deatoolbox.com
    
    #   Version: 1.0
#   LAST UPDATE: 26, April, 2017
    
    # Check size
    if size(X,1) != size(Y,1):
        error('Number of rows in X must be equal to number of rows in Y')
    
    
    # Get number of DMUs (n), inputs (m) and outputs (s)
    n,m=size(X,nargout=2)
# .\deaboot.m:39
    s=size(Y,2)
# .\deaboot.m:40
    
    options=getDEAoptions(n,varargin[arange()])
# .\deaboot.m:43
    
    if cellarray(['io','input']) == (options.orient):
        orient='io'
# .\deaboot.m:48
    else:
        if cellarray(['oo','output']) == (options.orient):
            orient='oo'
# .\deaboot.m:50
        else:
            if cellarray(['ddf']) == (options.orient):
                # orient = 'ddf';
                error('DEA Bootstrap not availalbe for \'ddf\'')
            else:
                if cellarray(['none']) == (options.orient):
                    error('Radial Model is Oriented')
    
    
    # RETURNS TO SCALE
    rts=options.rts
# .\deaboot.m:59
    
    if logical_not(isempty(options.Xeval)) and size(options.Xeval) != size(X):
        error('Xeval and X must be equal')
    
    
    if logical_not(isempty(options.Yeval)) and size(options.Yeval) != size(Y):
        error('Yeval and Y must be equal')
    
    
    # Original efficiency estimates
    if isempty(options.effRef):
        efforig=dea(X,Y,varargin[arange()])
# .\deaboot.m:72
        efforig=efforig.eff
# .\deaboot.m:73
    else:
        # Use reference efficiencies (when computing the RTS test)
        efforig=options.effRef
# .\deaboot.m:76
    
    # Invert efficiencies if 'io'
    if strcmp(orient,'io'):
        efforig_this=1 / efforig
# .\deaboot.m:81
    else:
        efforig_this=copy(efforig)
# .\deaboot.m:83
    
    # Get only inefficient units
    e=1e-05
# .\deaboot.m:87
    effn=efforig(efforig_this > 1 + e)
# .\deaboot.m:88
    
    eff2m=concat([[2 - effn],[effn]])
# .\deaboot.m:91
    
    #hm = 1.06 * min([eff2m_s, eff2m_iqr ./ 1.34]) * (length(eff2m)) ^(-1/5);
    hm=dot(dot(0.9,min(concat([std(eff2m),iqr(eff2m) / 1.34]))),(length(eff2m)) ** (- 1 / 5))
# .\deaboot.m:96
    
    # h = hm .*(length(eff2m) ./ length(effn)) * (std(efforig_this) ./ eff2m_s )
    h=dot(multiply(hm,((length(eff2m) / length(efforig))) ** (1 / 5)),(std(efforig_this) / std(eff2m)))
# .\deaboot.m:100
    
    nreps=options.nreps
# .\deaboot.m:103
    alph=options.alpha
# .\deaboot.m:104
    Eflag=nan(n,nreps)
# .\deaboot.m:106
    effBoot=nan(n,nreps)
# .\deaboot.m:107
    
    for i in arange(1,nreps).reshape(-1):
        # Dario and Simar (2007)
        # Invert efficiencies if 'io'
        efforigb=copy(efforig_this)
# .\deaboot.m:114
        eff2m=concat([[2 - efforigb],[efforigb]])
# .\deaboot.m:117
        effn=datasample(eff2m,n)
# .\deaboot.m:118
        effn_star=effn + multiply(h,randn(n,1))
# .\deaboot.m:121
        # eff_star = mean(effn) + (effn_star - mean(effn)) ./ (sqrt(1 + h^2 ./ var(effn)));
        eff_star=mean(effn) + (effn_star - mean(effn)) / (sqrt(1 + h ** 2 / var(eff2m)))
# .\deaboot.m:125
        eff_starr=copy(eff_star)
# .\deaboot.m:128
        eff_starr[eff_star < 1]=2 - eff_star(eff_star < 1)
# .\deaboot.m:129
        if strcmp(orient,'io'):
            eff_starr=1 / eff_starr
# .\deaboot.m:133
        # [5] Generate inefficient inputs or output and perform DEA
        Xref=[]
# .\deaboot.m:137
        Yref=[]
# .\deaboot.m:138
        if 'io' == (orient):
            Xref=multiply(repmat(efforig / eff_starr,1,m),X)
# .\deaboot.m:141
            Yref=copy(Y)
# .\deaboot.m:142
        else:
            if 'oo' == (orient):
                Xref=copy(X)
# .\deaboot.m:144
                Yref=multiply(repmat(efforig / eff_starr,1,s),Y)
# .\deaboot.m:145
        # Perform DEA with the subsample
        b=dea(Xref,Yref,'orient',orient,varargin[arange()],'Xeval',X,'Yeval',Y,'secondstep',0)
# .\deaboot.m:149
        effBoot[arange(),i]=b.eff
# .\deaboot.m:152
        Eflag[arange(),i]=b.exitflag(arange(),1)
# .\deaboot.m:154
    
    # Invert efficiencies if 'io'
    if strcmp(orient,'io'):
        effBoot=1 / effBoot
# .\deaboot.m:160
    
    # Bootstrap efficiency
    eff.bias = copy(mean(effBoot,2) - efforig_this)
# .\deaboot.m:164
    eff.b = copy(efforig_this - eff.bias)
# .\deaboot.m:165
    
    eff.c = copy(repelem(efforig_this,1,2) + quantile(repmat(efforig_this,1,nreps) - effBoot,concat([dot(0.5,alph),1 - dot(0.5,alph)]),2))
# .\deaboot.m:168
    
    eff.o = copy(efforig)
# .\deaboot.m:171
    eff.Boot = copy(effBoot)
# .\deaboot.m:172
    
    if strcmp(orient,'io'):
        eff.Boot = copy(1 / effBoot)
# .\deaboot.m:176
        eff.b = copy(1 / eff.b)
# .\deaboot.m:177
        eff.bias = copy(efforig - eff.b)
# .\deaboot.m:178
        eff.c = copy(1.0 / eff.c)
# .\deaboot.m:179
        eff.c = copy(eff.c(arange(),concat([2,1])))
# .\deaboot.m:180
    
    # Coompute variance
    eff.var = copy(var(effBoot,0,2))
# .\deaboot.m:184
    
    neval=copy(NaN)
# .\deaboot.m:187
    lambda_=copy(NaN)
# .\deaboot.m:188
    slack.X = copy(NaN)
# .\deaboot.m:189
    slack.Y = copy(NaN)
# .\deaboot.m:190
    Xeff=copy(NaN)
# .\deaboot.m:191
    Yeff=copy(NaN)
# .\deaboot.m:192
    
    out=deaout('n',n,'neval',neval.T,'s',s,'m',m,'X',X,'Y',Y,'names',options.names,'model','radial-bootstrap','orient',orient,'rts',rts,'lambda',lambda_,'slack',slack,'eff',eff,'Xeff',Xeff,'Yeff',Yeff,'exitflag',Eflag,'dispstr','names/eff.o/eff.b/eff.c','nreps',nreps,'alpha',alph)
# .\deaboot.m:195
    return out
    
if __name__ == '__main__':
    pass
    