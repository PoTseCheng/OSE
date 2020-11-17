# Generated with SMOP  0.41
from libsmop import *
# .\dea2table.m

    
@function
def dea2table(out=None,dispstr=None,*args,**kwargs):
    varargin = dea2table.varargin
    nargin = dea2table.nargin

    #DEA2TABLE Convert 'deaout' results into a table object
#   Convert data envelopment analysis results in a 'deaout' structure into
#   a MATLAB table object.
    
    #   T = DEA2TABLE( out ) Converts 'deaout' structure into a MATLAB table
#   object.
#   T = DEA2TABLE( out, dispstr ) Converts 'deaout' structure into a MATLAB
#   table object using the specified 'dispstr' structure.
    
    #   Example
#       
#       io = dea(X, Y, 'orient', 'io');
#       T = dea2table(io);
    
    #       T2 = dea2table(io, 'names/lambda/eff');
    
    #   See also DEAOUT, DEADISP
    
    #   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
#   http://www.deatoolbox.com
    
    #   Version: 1.0
#   LAST UPDATE: 9, May, 2017
    
    if nargin < 2:
        dispstr=out.dispstr
# .\dea2table.m:28
    
    
    # Convert to table
    dispstr=strsplit(dispstr,'/')
# .\dea2table.m:32
    
    T=table()
# .\dea2table.m:35
    
    for i in arange(1,length(dispstr)).reshape(-1):
        # Get param name
        paramstr=char(dispstr(i))
# .\dea2table.m:40
        dat=eval(sprintf('out.%s',paramstr))
# .\dea2table.m:43
        T=concat([T,table(dat)])
# .\dea2table.m:46
        name,__=getDEAformat(paramstr,out.orient,nargout=2)
# .\dea2table.m:49
        if isempty(name):
            disptext_field=sprintf('disptext_%s',strrep(paramstr,'.','_'))
# .\dea2table.m:53
            if isfield(out,disptext_field):
                # If custom name exists in the output structure use it
                name=eval(sprintf('out.%s',disptext_field))
# .\dea2table.m:56
            else:
                # If not, display paramstr name without eff.
                name=strrep(paramstr,'eff.','')
# .\dea2table.m:59
        # Store variable name
        T.Properties.VariableNames[size(T,2)]=cellstr(name)
# .\dea2table.m:64
    
    
    return T
    
if __name__ == '__main__':
    pass
    