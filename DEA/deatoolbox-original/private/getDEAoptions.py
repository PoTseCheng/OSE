# Generated with SMOP  0.41
from libsmop import *
# .\getDEAoptions.m

    
@function
def getDEAoptions(n=None,varargin=None,*args,**kwargs):
    varargin = getDEAoptions.varargin
    nargin = getDEAoptions.nargin

    #GETDEAOPTIONS Private function
#   Private function
    
    #   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
#   http://www.deatoolbox.com
    
    #   Version: 1.0
#   LAST UPDATE: 6, May, 2017
    
    
    # Default optimization-options
    optimoptsdef=optimoptions('linprog','display','off','Algorithm','dual-simplex','TolFun',1e-10,'TolCon',1e-07)
# .\getDEAoptions.m:13
    
    p=copy(inputParser)
# .\getDEAoptions.m:16
    if verLessThan('matlab','8.2'):
        addPar=lambda v1=None,v2=None,v3=None,v4=None: addParamValue(v1,v2,v3,v4)
# .\getDEAoptions.m:18
    else:
        addPar=lambda v1=None,v2=None,v3=None,v4=None: addParameter(v1,v2,v3,v4)
# .\getDEAoptions.m:20
    
    
    # Generic options
    addPar(p,'names',cellstr(int2str((arange(1,n)).T)),lambda x=None: iscellstr(x) and (length(x) == n))
    addPar(p,'optimopts',optimoptsdef,lambda x=None: logical_not(isempty(x)))
    addPar(p,'disp',0,lambda x=None: ismember(x,concat([0,1])))
    
    addPar(p,'orient','none',lambda x=None: any(validatestring(x,cellarray(['io','oo','ddf','none','ddf_cfg']))))
    addPar(p,'rts','crs',lambda x=None: any(validatestring(x,cellarray(['crs','vrs']))))
    addPar(p,'Gx',[],lambda x=None: isnumeric(x))
    addPar(p,'Gy',[],lambda x=None: isnumeric(x))
    addPar(p,'secondstep',1,lambda x=None: isnumeric(x))
    
    addPar(p,'Xeval',[],lambda x=None: isnumeric(x))
    addPar(p,'Yeval',[],lambda x=None: isnumeric(x))
    
    addPar(p,'rhoX',[],lambda x=None: isnumeric(x))
    addPar(p,'rhoY',[],lambda x=None: isnumeric(x))
    
    addPar(p,'fixbaset',[],lambda x=None: ismember(x,concat([0,1])))
    addPar(p,'geomean',[],lambda x=None: ismember(x,concat([0,1])))
    addPar(p,'period','geomean',lambda x=None: any(validatestring(x,cellarray(['base','comparison','geomean']))))
    
    addPar(p,'Xprice',[],lambda x=None: isnumeric(x))
    addPar(p,'Yprice',[],lambda x=None: isnumeric(x))
    
    addPar(p,'Yueval',[],lambda x=None: isnumeric(x))
    
    addPar(p,'nreps',200,lambda x=None: logical_and(isnumeric(x),(x > 0)))
    addPar(p,'alpha',0.05,lambda x=None: logical_and(isnumeric(x),(x > 0)))
    addPar(p,'effRef',[],lambda x=None: isnumeric(x))
    
    addPar(p,'warning',1,lambda x=None: ismember(x,concat([0,1])))
    p.parse(varargin[arange()])
    options=p.Results
# .\getDEAoptions.m:60
    
    if size(options.names,2) > 1:
        options.names = copy(options.names.T)
# .\getDEAoptions.m:64
    
    
    return options
    
if __name__ == '__main__':
    pass
    