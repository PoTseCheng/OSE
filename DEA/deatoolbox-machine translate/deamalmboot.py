# Generated with SMOP  0.41
from libsmop import *
# .\deamalmboot.m

    
@function
def deamalmboot(X=None,Y=None,varargin=None,*args,**kwargs):
    varargin = deamalmboot.varargin
    nargin = deamalmboot.nargin

    #DEAMALMBOOT Data envelopment analysis Malmquist indices bootstrap
#   Computes data envelopment analysis Malmquist indices bootstrap
    
    #   out = DEAMALMBOOT(X, Y, Name, Value) computes data envelopment analysis 
#   Malmquist indices bootstrap with inputs X and outputs Y. Model properties
#   are specified using one or more Name ,Value pair arguments.
    
    #   Additional properties:
#   - 'orient': orientation. Input oriented 'io', output oriented 'oo'.
#   - 'names': DMU names.
#   - 'fixbaset': previous year 0 (default), first year 1.
#   - 'nreps': number of bootstrap replications. Default is 200.
#   - 'alpha': alpha value for confidence intervals. Default is 0.05.
#   - 'period': compute geometric mean of base and comparison periods for 
#     technological change ('geomean'), use base period as reference ('base'),
#     or use comparison period as reference ('comparison').
#
    
    #   Deprecated parameters:
#   - 'geomean': compute geometric mean for technological change. Default
#     is 1. 'geomean' parameter has been deprecated and will dissapear in a
#     future realse. Set the new 'period' parapeter to 'geomean' for the 
#     previous behavior of 'geomean' = 1. Set 'period' to 'base' for the 
#     preivous behaviour of 'geomean' = 0.
    
    #   Example
#     
#      iomalm = deamalmboot(X, Y, 'orient', 'io', 'nreps, 200);
    
    #   See also DEAOUT, DEA, DEABOOT, DEAMALML
    
    #   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
#   http://www.deatoolbox.com
    
    #   Version: 1.0
#   LAST UPDATE: 11, July, 2018
    
    # Check size
    if size(X,1) != size(Y,1):
        error('Number of rows in X must be equal to number of rows in Y')
    
    
    if size(X,3) != size(Y,3):
        error('Number of time periods in X and Y must be equal')
    
    
    # Get number of DMUs (n), inputs (m) and outputs (s)
    n,m,T=size(X,nargout=3)
# .\deamalmboot.m:50
    s=size(Y,2)
# .\deamalmboot.m:51
    
    options=getDEAoptions(n,varargin[arange()])
# .\deamalmboot.m:54
    
    if logical_not(isempty(options.Xeval)) and size(options.Xeval) != size(X):
        error('Xeval and X must be equal')
    
    
    if logical_not(isempty(options.Yeval)) and size(options.Yeval) != size(Y):
        error('Yeval and Y must be equal')
    
    
    # Check RTS
    if logical_not(strcmp(options.rts,'crs')):
        error('Malmquist index only available for \'crs\' returns to scale')
    
    
    # Check orientation
    if strcmp(options.orient,'ddf') or strcmp(options.orient,'oo'):
        error('Malmquist index bootstrap only for \'io\' orientation')
    
    
    # Get number of Bootstrap replications and significance
    nreps=options.nreps
# .\deamalmboot.m:76
    alph=options.alpha
# .\deamalmboot.m:77
    
    Mb=nan(n,T - 1)
# .\deamalmboot.m:80
    MTECb=nan(n,T - 1)
# .\deamalmboot.m:81
    MTCb=nan(n,T - 1)
# .\deamalmboot.m:82
    MB=nan(n,nreps,T - 1)
# .\deamalmboot.m:84
    MTECB=nan(n,nreps,T - 1)
# .\deamalmboot.m:85
    MTCB=nan(n,nreps,T - 1)
# .\deamalmboot.m:86
    EflagB=nan(n,dot(nreps,2),(T - 1))
# .\deamalmboot.m:87
    
    
    # Check if 'geomean' and the old parameter 'period' are correct
    if logical_not(isempty(options.geomean)):
        warning('\'geomean\' parameter has been deprecated and will dissapear in a future realse.\\n Set the new \'period\' parameter to \'geomean\' for the previous behavior of \'geomean\' = 1.\\n Set \'period\' to \'base\' for the preivous behaviour of \'geomean\' = 0. See help for more information.','DEATOOLBOX:deprecated')
        if options.geomean:
            if logical_not(strcmp(options.period,'geomean')):
                error('If \'geomean\' is set to 1, \'period\' must be set to \'geomean\'')
        else:
            if logical_not(strcmp(options.period,'base')):
                error('If \'geomean\' is set to 0, \'period\' must be set to \'base\'')
    
    
    # Original Malmquist indices
    tempmalm=deamalm(X,Y,varargin[arange()])
# .\deamalmboot.m:104
    Mo=tempmalm.eff.M
# .\deamalmboot.m:105
    MTECo=tempmalm.eff.MTEC
# .\deamalmboot.m:106
    MTCo=tempmalm.eff.MTC
# .\deamalmboot.m:107
    
    # Silverman's (1986) suggestion for bivariate data
    h=(4 / (multiply(5,n))) ** (1 / 6)
# .\deamalmboot.m:111
    
    for t in arange(1,T - 1).reshape(-1):
        # Get base period
        if isempty(options.fixbaset) or options.fixbaset == 0:
            tb=copy(t)
# .\deamalmboot.m:117
        else:
            if options.fixbaset == 1:
                tb=1
# .\deamalmboot.m:119
        # A, B, and Delta matrix
        A=dea(X(arange(),arange(),tb),Y(arange(),arange(),tb),varargin[arange()],'secondstep',0)
# .\deamalmboot.m:123
        A=A.eff
# .\deamalmboot.m:124
        B=dea(X(arange(),arange(),t + 1),Y(arange(),arange(),t + 1),varargin[arange()],'secondstep',0)
# .\deamalmboot.m:125
        B=B.eff
# .\deamalmboot.m:126
        Delta=concat([[A,B],[2 - A,B],[2 - A,2 - B],[A,2 - B]])
# .\deamalmboot.m:127
        DeltaStarReps=nan(n,2,nreps)
# .\deamalmboot.m:133
        idxReps=nan(1,n,nreps)
# .\deamalmboot.m:134
        for i in arange(1,nreps).reshape(-1):
            DeltaStarReps(arange(),arange(),i),idxReps(arange(),arange(),i)=datasample(Delta,n,nargout=2)
# .\deamalmboot.m:136
        # For each replication
        for i in arange(1,nreps).reshape(-1):
            # Get random sample
            DeltaStar=DeltaStarReps(arange(),arange(),i)
# .\deamalmboot.m:143
            idx=idxReps(arange(),arange(),i)
# .\deamalmboot.m:144
            deltaMat=diag(mean(DeltaStar))
# .\deamalmboot.m:147
            Sigma=cov(A,B)
# .\deamalmboot.m:150
            SigmaR=copy(Sigma)
# .\deamalmboot.m:151
            SigmaR[1,2]=- SigmaR(1,2)
# .\deamalmboot.m:152
            SigmaR[2,1]=- SigmaR(2,1)
# .\deamalmboot.m:153
            raMat=nan(n,2)
# .\deamalmboot.m:156
            idxSigma=idx <= logical_or(n,(idx > logical_and((dot(2,n)),idx) <= dot(3,n)))
# .\deamalmboot.m:158
            raMat[idxSigma,arange()]=mvnrnd(concat([0,0]),Sigma,sum(idxSigma))
# .\deamalmboot.m:159
            idxSigmaR=logical_or((idx > logical_and(n,idx) <= dot(2,n)),idx) > dot(3,n)
# .\deamalmboot.m:161
            raMat[idxSigmaR,arange()]=mvnrnd(concat([0,0]),SigmaR,sum(idxSigmaR))
# .\deamalmboot.m:162
            C=ones(n,2)
# .\deamalmboot.m:166
            Gamma=multiply((1 + h ** 2) ** (- 1 / 2),(DeltaStar + multiply(h,raMat) - dot(C,deltaMat))) + (dot(C,deltaMat))
# .\deamalmboot.m:167
            GammaStar=copy(Gamma)
# .\deamalmboot.m:172
            GammaStar[Gamma < 1]=2 - Gamma(Gamma < 1)
# .\deamalmboot.m:173
            Xpseudo1=multiply(X(arange(),arange(),tb),repmat(GammaStar(arange(),1) / A,1,m))
# .\deamalmboot.m:176
            Xpseudo2=multiply(X(arange(),arange(),t + 1),repmat(GammaStar(arange(),2) / B,1,m))
# .\deamalmboot.m:177
            # Compute efficiency at base period
            temp_dea=dea(Xpseudo1,Y(arange(),arange(),tb),varargin[arange()],'secondstep',0,'Xeval',X(arange(),arange(),tb),'Yeval',Y(arange(),arange(),tb))
# .\deamalmboot.m:181
            tb_eff=temp_dea.eff
# .\deamalmboot.m:182
            temp_dea=dea(Xpseudo2,Y(arange(),arange(),t + 1),varargin[arange()],'secondstep',0,'Xeval',X(arange(),arange(),t + 1),'Yeval',Y(arange(),arange(),t + 1))
# .\deamalmboot.m:185
            t1_eff=temp_dea.eff
# .\deamalmboot.m:186
            temp_dea=dea(Xpseudo1,Y(arange(),arange(),tb),varargin[arange()],'Xeval',Xpseudo2,'Yeval',Y(arange(),arange(),t + 1),'secondstep',0)
# .\deamalmboot.m:189
            tbevalt1_eff=temp_dea.eff
# .\deamalmboot.m:194
            t1evaltb_eff=copy(NaN)
# .\deamalmboot.m:197
            if cellarray(['geomean','comparison']) == (options.period):
                # Evaluate each DMU at t + 1, with the others at base period
                temp_dea=dea(Xpseudo2,Y(arange(),arange(),t + 1),varargin[arange()],'Xeval',Xpseudo1,'Yeval',Y(arange(),arange(),tb),'secondstep',0)
# .\deamalmboot.m:201
                t1evaltb_eff=temp_dea.eff
# .\deamalmboot.m:205
            # Technical Efficiency
            MTECB[arange(),i,t]=t1_eff / tb_eff
# .\deamalmboot.m:209
            if 'geomean' == (options.period):
                MTCB[arange(),i,t]=(multiply((tbevalt1_eff / t1_eff),(tb_eff / t1evaltb_eff))) ** (1 / 2)
# .\deamalmboot.m:214
            else:
                if 'base' == (options.period):
                    MTCB[arange(),i,t]=tbevalt1_eff / t1_eff
# .\deamalmboot.m:216
                else:
                    if 'comparison' == (options.period):
                        MTCB[arange(),i,t]=tb_eff / t1evaltb_eff
# .\deamalmboot.m:218
            # Malmquist index
            MB[arange(),i,t]=multiply(MTECB(arange(),i,t),MTCB(arange(),i,t))
# .\deamalmboot.m:222
        # Bootrstrap Technical Efficiency
        MTEC.bias[arange(),t]=mean(MTECB(arange(),arange(),t),2) - MTECo(arange(),t)
# .\deamalmboot.m:227
        MTEC.b[arange(),t]=MTECo(arange(),t) - MTEC.bias(arange(),t)
# .\deamalmboot.m:228
        confInt=repelem(MTECo(arange(),t),1,2) + quantile(repmat(MTECo(arange(),t),1,nreps) - MTECB(arange(),arange(),t),concat([dot(0.5,alph),1 - dot(0.5,alph)]),2)
# .\deamalmboot.m:229
        MTEC.cL[arange(),t]=confInt(arange(),1)
# .\deamalmboot.m:231
        MTEC.cU[arange(),t]=confInt(arange(),2)
# .\deamalmboot.m:232
        MTC.bias[arange(),t]=mean(MTCB(arange(),arange(),t),2) - MTCo(arange(),t)
# .\deamalmboot.m:235
        MTC.b[arange(),t]=MTCo(arange(),t) - MTC.bias(arange(),t)
# .\deamalmboot.m:236
        confInt=repelem(MTCo(arange(),t),1,2) + quantile(repmat(MTCo(arange(),t),1,nreps) - MTCB(arange(),arange(),t),concat([dot(0.5,alph),1 - dot(0.5,alph)]),2)
# .\deamalmboot.m:237
        MTC.cL[arange(),t]=confInt(arange(),1)
# .\deamalmboot.m:239
        MTC.cU[arange(),t]=confInt(arange(),2)
# .\deamalmboot.m:240
        M.bias[arange(),t]=mean(MB(arange(),arange(),t),2) - Mo(arange(),t)
# .\deamalmboot.m:243
        M.b[arange(),t]=Mo(arange(),t) - M.bias(arange(),t)
# .\deamalmboot.m:244
        confInt=repelem(Mo(arange(),t),1,2) + quantile(repmat(Mo(arange(),t),1,nreps) - MB(arange(),arange(),t),concat([dot(0.5,alph),1 - dot(0.5,alph)]),2)
# .\deamalmboot.m:245
        M.cL[arange(),t]=confInt(arange(),1)
# .\deamalmboot.m:247
        M.cU[arange(),t]=confInt(arange(),2)
# .\deamalmboot.m:248
    
    
    # Store original malmquist
    M.o = copy(Mo)
# .\deamalmboot.m:253
    MTEC.o = copy(MTECo)
# .\deamalmboot.m:254
    MTC.o = copy(MTCo)
# .\deamalmboot.m:255
    
    eff.M = copy(M)
# .\deamalmboot.m:258
    eff.MTEC = copy(MTEC)
# .\deamalmboot.m:259
    eff.MTC = copy(MTC)
# .\deamalmboot.m:260
    eff.T = copy(T)
# .\deamalmboot.m:261
    
    neval=copy(NaN)
# .\deamalmboot.m:264
    lambda_=copy(NaN)
# .\deamalmboot.m:265
    slack.X = copy(NaN)
# .\deamalmboot.m:266
    slack.Y = copy(NaN)
# .\deamalmboot.m:267
    Xeff=copy(NaN)
# .\deamalmboot.m:268
    Yeff=copy(NaN)
# .\deamalmboot.m:269
    Eflag=copy(NaN)
# .\deamalmboot.m:270
    
    out=deaout('n',n,'neval',neval.T,'s',s,'m',m,'X',X,'Y',Y,'names',options.names,'model','radial-malmquist-bootstrap','orient',options.orient,'rts',options.rts,'lambda',lambda_,'slack',slack,'eff',eff,'Xeff',Xeff,'Yeff',Yeff,'exitflag',Eflag,'dispstr','names/eff.M.o/eff.M.b/eff.M.cL/eff.M.cU','nreps',nreps,'alpha',alph)
# .\deamalmboot.m:273
    out.period = copy(options.period)
# .\deamalmboot.m:282
    out.fixbaset = copy(options.fixbaset)
# .\deamalmboot.m:283
    
    out.disptext_text2 = copy('Malmquist:')
# .\deamalmboot.m:286
    out.disptext_text4 = copy('M = Malmquist. Mboot = Bootstrapped Malmquist. McLow = Lower confidence interval. McUpp: Upper confidence interval.')
# .\deamalmboot.m:287
    return out
    
if __name__ == '__main__':
    pass
    