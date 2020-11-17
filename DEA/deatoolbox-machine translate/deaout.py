# Generated with SMOP  0.41
from libsmop import *
# .\deaout.m

    
@function
def deaout(varargin=None,*args,**kwargs):
    varargin = deaout.varargin
    nargin = deaout.nargin

    #DEAOUT Generates a default deaout structure to store dea results
#   Generates a default deaout structure to store dea results
    
    #   Information fields:
#   - n: number of DMU's.
#   - neval: number of evaluated DMU's.
#   - m: number of inputs.
#   - s: number of outputs.
#   - r: number of undesirable outputs.
#   - model: dea model.
#   - orient: orientation.
#   - rts: returns to scale.
    
    #   Common fields:
#   - names: DMU names.
#   - X: inputs.
#   - Y: outputs.
#   - eff: efficiency measure
#   - slack.X: input slacks.
#   - slack.Y: ouput slacks.
#   - lambda: computed lambda'.
#   - Xeff: efficient X's.
#   - Yeff: efficient Y's.
#   - dual.X: input shadow prices.
#   - dual.Y: output shadow prices.
#   - dual.rts: RTS dual.
#   - exitflag: exit flags of the optimization.
    
    #   Scale efficiency models:
#   - eff.crs: CRS efficiency.
#   - eff.vrs: VRS efficiency.
#   - eff.scale: Scale efficiency.
    
    #   Malmquist index:
#   - eff.M: Malmquist index.
#   - eff.MTEC: Technical efficiency change.
#   - eff.MTC: Technical change.
    
    #   Allocative efficiency model:
#   - W: Cost.
#   - P: Price.
#   - eff.C: Cost efficiency.
#   - eff.R: Revenue efficiency.
#   - eff.P: Profit efficiency.
#   - eff.A: Allocative efficiency.
#   - eff.T: Technical efficiency.
    
    #   Undesirable outputs model;
#   - Yu: Undesirable outputs.
#   - slack.Yu: Undesirable outputs slacks.
    
    #   Malmquist-Luenberger index:
#   - eff.ML: Malmquist-Luenberger index.
#   - eff.MLTEC: Technical efficiency change.
#   - eff.MLTC: Technical change.
    
    #   DEA Bootstrap:
#   - eff.o: Origianl efficiency.
#   - eff.b: Bootstraped efficiency.
#   - eff.c: Efficiency confidence interval.
    
    #   Example
#     
#      model = deaout();
    
    #   See also DEA
    
    #   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
#   http://www.deatoolbox.com
    
    #   Version: 1.0
#   LAST UPDATE: 1, September, 2016
#
    
    
    # Parse options
    p=copy(inputParser)
# .\deaout.m:78
    if verLessThan('matlab','8.2'):
        addPar=lambda v1=None,v2=None,v3=None,v4=None: addParamValue(v1,v2,v3,v4)
# .\deaout.m:80
    else:
        addPar=lambda v1=None,v2=None,v3=None,v4=None: addParameter(v1,v2,v3,v4)
# .\deaout.m:82
    
    # Generic options
    addPar(p,'n',NaN,lambda x=None: isnumeric(x))
    addPar(p,'neval',NaN,lambda x=None: isnumeric(x))
    addPar(p,'s',NaN,lambda x=None: isnumeric(x))
    addPar(p,'m',NaN,lambda x=None: isnumeric(x))
    addPar(p,'X',NaN,lambda x=None: isnumeric(x))
    addPar(p,'Y',NaN,lambda x=None: isnumeric(x))
    addPar(p,'names',[],lambda x=None: iscellstr(x))
    addPar(p,'model','radial',lambda x=None: any(validatestring(x,cellarray(['radial','radial-supereff','radial-malmquist','directional','directional-supereff','directional-undesirable','directional-malmquist-luenberger','additive','additive-supereff','additive-profit','allocative-cost','allocative-revenue','allocative-profit','radial-bootstrap','radial-malmquist-bootstrap']))))
    addPar(p,'orient','none',lambda x=None: any(validatestring(x,cellarray(['io','oo','ddf','none','ddf_cfg']))))
    addPar(p,'rts','crs',lambda x=None: any(validatestring(x,cellarray(['crs','vrs','scaleeff']))))
    addPar(p,'lambda',NaN,lambda x=None: isnumeric(x))
    addPar(p,'slack',NaN,lambda x=None: isstruct(x))
    addPar(p,'eff',NaN,lambda x=None: isnumeric(x) or isstruct(x))
    addPar(p,'Xeff',NaN,lambda x=None: isnumeric(x))
    addPar(p,'Yeff',NaN,lambda x=None: isnumeric(x))
    addPar(p,'dual',NaN,lambda x=None: isstruct(x))
    addPar(p,'dispstr','names/X/Y/slackX/slackY/eff',lambda x=None: ischar(x))
    addPar(p,'exitflag',NaN,lambda x=None: isnumeric(x))
    
    addPar(p,'Xprice',NaN,lambda x=None: isnumeric(x))
    addPar(p,'Yprice',NaN,lambda x=None: isnumeric(x))
    
    addPar(p,'r',NaN,lambda x=None: isnumeric(x))
    addPar(p,'Yu',NaN,lambda x=None: isnumeric(x))
    addPar(p,'Yueff',NaN,lambda x=None: isnumeric(x))
    
    addPar(p,'nreps',NaN,lambda x=None: isnumeric(x))
    addPar(p,'alpha',NaN,lambda x=None: isnumeric(x))
    p.parse(varargin[arange()])
    out=p.Results
# .\deaout.m:126
    return out
    
if __name__ == '__main__':
    pass
    