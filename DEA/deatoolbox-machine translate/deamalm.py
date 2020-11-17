# Generated with SMOP  0.41
from libsmop import *
# .\deamalm.m

    
@function
def deamalm(X=None,Y=None,varargin=None,*args,**kwargs):
    varargin = deamalm.varargin
    nargin = deamalm.nargin

    #DEAMALM Data envelopment analysis Malmquist indices
#   Computes data envelopment analysis Malmquist indices
    
    #   out = DEAMALM(X, Y, Name, Value) computes data envelopment analysis 
#   Malmquist indices with inputs X and outputs Y. Model properties are 
#   specified using one or more Name ,Value pair arguments.
    
    #   Additional properties:
#   - 'orient': orientation. Input oriented 'io', output oriented 'oo'.
#   - 'names': DMU names.
#   - 'fixbaset': base year is previous year 0 (default), or always the 
#     first year 1.
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
#      iomalm = deamalm(X, Y, 'orient', 'io');
#      oomalmfixex = deamalm(X, Y, 'orient', 'oo', 'fixbaset', 1);
    
    #   See also DEAOUT, DEA, DEAMALMLUEN
    
    #   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
#   http://www.deatoolbox.com
    
    #   Version: 1.0
#   LAST UPDATE: 6, May, 2017
    
    # Check size
    if size(X,1) != size(Y,1):
        error('Number of rows in X must be equal to number of rows in Y')
    
    
    if size(X,3) != size(Y,3):
        error('Number of time periods in X and Y must be equal')
    
    
    # Get number of DMUs (n), inputs (m) and outputs (s)
    n,m,T=size(X,nargout=3)
# .\deamalm.m:50
    s=size(Y,2)
# .\deamalm.m:51
    
    options=getDEAoptions(n,varargin[arange()])
# .\deamalm.m:54
    
    if logical_not(isempty(options.Xeval)) and size(options.Xeval) != size(X):
        error('Xeval and X must be equal')
    
    
    if logical_not(isempty(options.Yeval)) and size(options.Yeval) != size(Y):
        error('Yeval and Y must be equal')
    
    
    # Check orientation
    if strcmp(options.orient,'ddf'):
        error('Malmquist index for \'ddf\' not yet implemented')
    
    
    # Create matrices to store results
    M=nan(n,T - 1)
# .\deamalm.m:71
    MTEC=nan(n,T - 1)
# .\deamalm.m:72
    MTC=nan(n,T - 1)
# .\deamalm.m:73
    Eflag=nan(n,dot((T - 1),4))
# .\deamalm.m:74
    
    if logical_not(isempty(options.geomean)):
        warning('\'geomean\' parameter has been deprecated and will dissapear in a future realse.\\n Set the new \'period\' parameter to \'geomean\' for the previous behavior of \'geomean\' = 1.\\n Set \'period\' to \'base\' for the preivous behaviour of \'geomean\' = 0. See help for more information.','DEATOOLBOX:deprecated')
        if options.geomean:
            if logical_not(strcmp(options.period,'geomean')):
                error('If \'geomean\' is set to 1, \'period\' must be set to \'geomean\'')
        else:
            if logical_not(strcmp(options.period,'base')):
                error('If \'geomean\' is set to 0, \'period\' must be set to \'base\'')
    
    
    # For each time period
    for t in arange(1,T - 1).reshape(-1):
        # Get base period
        if isempty(options.fixbaset) or options.fixbaset == 0:
            tb=copy(t)
# .\deamalm.m:94
        else:
            if options.fixbaset == 1:
                tb=1
# .\deamalm.m:96
        # Compute efficiency at base period
        temp_dea=dea(X(arange(),arange(),tb),Y(arange(),arange(),tb),varargin[arange()],'secondstep',0)
# .\deamalm.m:100
        tb_eff=temp_dea.eff
# .\deamalm.m:101
        Eflag[arange(),(dot(2,t) - 1)]=temp_dea.exitflag(arange(),1)
# .\deamalm.m:102
        temp_dea=dea(X(arange(),arange(),t + 1),Y(arange(),arange(),t + 1),varargin[arange()],'secondstep',0)
# .\deamalm.m:105
        t1_eff=temp_dea.eff
# .\deamalm.m:106
        Eflag[arange(),(dot(2,t) - 1) + 1]=temp_dea.exitflag(arange(),1)
# .\deamalm.m:107
        temp_dea=dea(X(arange(),arange(),tb),Y(arange(),arange(),tb),varargin[arange()],'Xeval',X(arange(),arange(),t + 1),'Yeval',Y(arange(),arange(),t + 1),'secondstep',0)
# .\deamalm.m:110
        tbevalt1_eff=temp_dea.eff
# .\deamalm.m:114
        Eflag[arange(),(dot(2,t) - 1) + 2]=temp_dea.exitflag(arange(),1)
# .\deamalm.m:115
        if cellarray(['geomean','comparison']) == (options.period):
            # Evaluate each DMU at base period, with the others at t + 1
            temp_dea=dea(X(arange(),arange(),t + 1),Y(arange(),arange(),t + 1),varargin[arange()],'Xeval',X(arange(),arange(),tb),'Yeval',Y(arange(),arange(),tb),'secondstep',0)
# .\deamalm.m:121
            t1evaltb_eff=temp_dea.eff
# .\deamalm.m:125
            Eflag[arange(),(dot(2,t) - 1) + 3]=temp_dea.exitflag(arange(),1)
# .\deamalm.m:126
        else:
            if 'base' == (options.period):
                t1evaltb_eff=copy(NaN)
# .\deamalm.m:128
        # Inverse efficiencies if 'oo'
        if strcmp(options.orient,'oo'):
            tb_eff=1 / tb_eff
# .\deamalm.m:133
            t1_eff=1 / t1_eff
# .\deamalm.m:134
            tbevalt1_eff=1 / tbevalt1_eff
# .\deamalm.m:135
            t1evaltb_eff=1 / t1evaltb_eff
# .\deamalm.m:136
        # Technical Efficiency
        MTEC[arange(),t]=t1_eff / tb_eff
# .\deamalm.m:140
        if 'geomean' == (options.period):
            MTC[arange(),t]=(multiply((tbevalt1_eff / t1_eff),(tb_eff / t1evaltb_eff))) ** (1 / 2)
# .\deamalm.m:145
        else:
            if 'base' == (options.period):
                MTC[arange(),t]=tbevalt1_eff / t1_eff
# .\deamalm.m:147
            else:
                if 'comparison' == (options.period):
                    MTC[arange(),t]=tb_eff / t1evaltb_eff
# .\deamalm.m:149
        # Malmquist index
        M[arange(),t]=multiply(MTEC(arange(),t),MTC(arange(),t))
# .\deamalm.m:153
    
    
    # Store Malmquist results in the efficiency structure
    eff.M = copy(M)
# .\deamalm.m:158
    eff.MTEC = copy(MTEC)
# .\deamalm.m:159
    eff.MTC = copy(MTC)
# .\deamalm.m:160
    eff.T = copy(T)
# .\deamalm.m:161
    
    neval=copy(NaN)
# .\deamalm.m:164
    lambda_=copy(NaN)
# .\deamalm.m:165
    slack.X = copy(NaN)
# .\deamalm.m:166
    slack.Y = copy(NaN)
# .\deamalm.m:167
    Xeff=copy(NaN)
# .\deamalm.m:168
    Yeff=copy(NaN)
# .\deamalm.m:169
    
    out=deaout('n',n,'neval',neval.T,'s',s,'m',m,'X',X,'Y',Y,'names',options.names,'model','radial-malmquist','orient',options.orient,'rts',options.rts,'lambda',lambda_,'slack',slack,'eff',eff,'Xeff',Xeff,'Yeff',Yeff,'exitflag',Eflag,'dispstr','names/eff.M/eff.MTEC/eff.MTC')
# .\deamalm.m:172
    out.period = copy(options.period)
# .\deamalm.m:180
    out.fixbaset = copy(options.fixbaset)
# .\deamalm.m:181
    
    out.disptext_text2 = copy('Malmquist:')
# .\deamalm.m:184
    out.disptext_text4 = copy('M = Malmquist. MTEC = Technical Efficiency Change. MTC = Technical Change.')
# .\deamalm.m:185
    return out
    
if __name__ == '__main__':
    pass
    