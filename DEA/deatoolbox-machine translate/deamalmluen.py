# Generated with SMOP  0.41
from libsmop import *
# .\deamalmluen.m

    
@function
def deamalmluen(X=None,Y=None,Yu=None,varargin=None,*args,**kwargs):
    varargin = deamalmluen.varargin
    nargin = deamalmluen.nargin

    #DEAMALMLUEN Data envelopment analysis Malmquist-Luenberger indices
#   Computes data envelopment analysis Malmquist-Luenberger indices
    
    #   out = DEAMALMLUEN(X, Y, Yu, Name, Value) computes data envelopment 
#   analysis Malmquist-Luenberger indices with inputs X and outputs Y. 
#   Model properties are specified using one or more Name ,Value pair 
#   arguments.
    
    #   Additional properties:
#   - 'names': DMU names.
#   - 'orient': orientation. Directional distane function with undesirable
#   outputs 'ddf' (Aparicio, Pastor and Zofio, 2013), default. Directional 
#   distance function with undesirable outputs 'ddf_ccf' (Chung, Fare and 
#   Grosskopf).
#   - 'fixbaset': previous year 0 (default), first year 1.
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
#      malmluen = deamalmluen(X, Y, Yu);
#      malmluenfixed = deamalmluen(X, Y, Yu, 'fixbaset', 1);
    
    #   See also DEAOUT, DEA, DEAUND, DEAMALMLUEN
    
    #   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
#   http://www.deatoolbox.com
    
    #   Version: 1.0
#   LAST UPDATE: 6, May, 2017
    
    # Check size
    if size(X,1) != size(Y,1):
        error('Number of rows in X must be equal to number of rows in Y')
    
    
    if size(X,3) != size(Y,3):
        error('Number of time periods in X and Y must be equal')
    
    
    if size(Yu,3) != size(Y,3):
        error('Number of time periods in Y and Yu must be equal')
    
    
    # Get number of DMUs (n), inputs (m), outputs (s), and undesirable
    # outputs (r)
    n,m,T=size(X,nargout=3)
# .\deamalmluen.m:58
    s=size(Y,2)
# .\deamalmluen.m:59
    r=size(Yu,2)
# .\deamalmluen.m:60
    
    options=getDEAoptions(n,varargin[arange()])
# .\deamalmluen.m:63
    
    if logical_not(isempty(options.Xeval)) and size(options.Xeval) != size(X):
        error('Xeval and X must be equal')
    
    
    if logical_not(isempty(options.Yeval)) and size(options.Yeval) != size(Y):
        error('Yeval and Y must be equal')
    
    
    # Check RTS
    if logical_not(strcmp(options.rts,'crs')):
        error('Malmquist-Luenberger index only available for \'crs\' returns to scale')
    
    
    # Replace default 'none' orientation to 'ddf'
    orient=options.orient
# .\deamalmluen.m:80
    if strcmp(orient,'none'):
        orient='ddf'
# .\deamalmluen.m:82
    
    
    # Check orientation
    if logical_not(strcmp(orient,'ddf')) and logical_not(strcmp(orient,'ddf_cfg')):
        error('Malmquist-Luenberger index is for \'ddf\' or \'ddf_ccf\' with undesirable outputs')
    
    
    # Create matrices to store results
    ML=nan(n,T - 1)
# .\deamalmluen.m:91
    MLTEC=nan(n,T - 1)
# .\deamalmluen.m:92
    MLTC=nan(n,T - 1)
# .\deamalmluen.m:93
    if options.geomean:
        Eflag=nan(n,dot((T - 1),4))
# .\deamalmluen.m:95
    else:
        Eflag=nan(n,dot((T - 1),3))
# .\deamalmluen.m:97
    
    
    # Check if 'geomean' and the old parameter 'period' are correct
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
# .\deamalmluen.m:118
        else:
            if options.fixbaset == 1:
                tb=1
# .\deamalmluen.m:120
        # Compute efficiency at base period
        temp_dea=deaund(X(arange(),arange(),tb),Y(arange(),arange(),tb),Yu(arange(),arange(),tb),varargin[arange()],'secondstep',0)
# .\deamalmluen.m:124
        tb_eff=temp_dea.eff
# .\deamalmluen.m:125
        Eflag[arange(),(dot(2,t) - 1)]=temp_dea.exitflag(arange(),1)
# .\deamalmluen.m:126
        temp_dea=deaund(X(arange(),arange(),t + 1),Y(arange(),arange(),t + 1),Yu(arange(),arange(),t + 1),varargin[arange()],'secondstep',0)
# .\deamalmluen.m:129
        t1_eff=temp_dea.eff
# .\deamalmluen.m:130
        Eflag[arange(),(dot(2,t) - 1) + 1]=temp_dea.exitflag(arange(),1)
# .\deamalmluen.m:131
        temp_dea=deaund(X(arange(),arange(),tb),Y(arange(),arange(),tb),Yu(arange(),arange(),tb),varargin[arange()],'Xeval',X(arange(),arange(),t + 1),'Yeval',Y(arange(),arange(),t + 1),'Yueval',Yu(arange(),arange(),t + 1),'secondstep',0)
# .\deamalmluen.m:134
        tbevalt1_eff=temp_dea.eff
# .\deamalmluen.m:138
        Eflag[arange(),(dot(2,t) - 1) + 2]=temp_dea.exitflag(arange(),1)
# .\deamalmluen.m:139
        if cellarray(['geomean','comparison']) == (options.period):
            # Evaluate each DMU at t + 1, with the others at base period
            temp_dea=deaund(X(arange(),arange(),t + 1),Y(arange(),arange(),t + 1),Yu(arange(),arange(),t + 1),varargin[arange()],'Xeval',X(arange(),arange(),tb),'Yeval',Y(arange(),arange(),tb),'Yueval',Yu(arange(),arange(),tb),'secondstep',0)
# .\deamalmluen.m:145
            t1evaltb_eff=temp_dea.eff
# .\deamalmluen.m:150
            Eflag[arange(),(dot(2,t) - 1) + 3]=temp_dea.exitflag(arange(),1)
# .\deamalmluen.m:151
        else:
            if 'base' == (options.period):
                t1evaltb_eff=copy(NaN)
# .\deamalmluen.m:153
        # Technical Efficiency
        MLTEC[arange(),t]=(1 + tb_eff) / (1 + t1_eff)
# .\deamalmluen.m:157
        if 'geomean' == (options.period):
            MLTC[arange(),t]=(multiply(((1 + t1_eff) / (1 + tbevalt1_eff)),((1 + t1evaltb_eff) / (1 + tb_eff)))) ** (1 / 2)
# .\deamalmluen.m:162
        else:
            if 'base' == (options.period):
                MLTC[arange(),t]=(1 + t1_eff) / (1 + tbevalt1_eff)
# .\deamalmluen.m:164
            else:
                if 'comparison' == (options.period):
                    MLTC[arange(),t]=(1 + t1evaltb_eff) / (1 + tb_eff)
# .\deamalmluen.m:166
        # Malmquist-Luenberger index
        ML[arange(),t]=multiply(MLTEC(arange(),t),MLTC(arange(),t))
# .\deamalmluen.m:170
    
    
    # Store Malmquist results in the efficiency structure
    eff.ML = copy(ML)
# .\deamalmluen.m:175
    eff.MLTEC = copy(MLTEC)
# .\deamalmluen.m:176
    eff.MLTC = copy(MLTC)
# .\deamalmluen.m:177
    eff.T = copy(T)
# .\deamalmluen.m:178
    
    neval=copy(NaN)
# .\deamalmluen.m:182
    lambda_=copy(NaN)
# .\deamalmluen.m:183
    slack.X = copy(NaN)
# .\deamalmluen.m:184
    slack.Y = copy(NaN)
# .\deamalmluen.m:185
    slack.Yu = copy(NaN)
# .\deamalmluen.m:186
    Xeff=copy(NaN)
# .\deamalmluen.m:187
    Yeff=copy(NaN)
# .\deamalmluen.m:188
    
    out=deaout('n',n,'neval',neval.T,'s',s,'m',m,'X',X,'Y',Y,'names',options.names,'model','directional-malmquist-luenberger','orient',orient,'rts',options.rts,'lambda',lambda_,'slack',slack,'eff',eff,'Xeff',Xeff,'Yeff',Yeff,'exitflag',Eflag,'dispstr','names/eff.ML/eff.MLTEC/eff.MLTC','r',r,'Yu',Yu)
# .\deamalmluen.m:191
    out.period = copy(options.period)
# .\deamalmluen.m:200
    out.fixbaset = copy(options.fixbaset)
# .\deamalmluen.m:201
    
    out.disptext_text2 = copy('Malmquist-Luenberger:')
# .\deamalmluen.m:204
    out.disptext_text4 = copy('ML: Malmquist-Luenberger. MLTEC: Technical Efficiency Change. MLTC: Technical Change.')
# .\deamalmluen.m:205
    return out
    
if __name__ == '__main__':
    pass
    