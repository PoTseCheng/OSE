# Generated with SMOP  0.41
from libsmop import *
# .\deadisp.m

    
@function
def deadisp(out=None,dispstr=None,*args,**kwargs):
    varargin = deadisp.varargin
    nargin = deadisp.nargin

    #DEADISP Display data envelopment analysis results
#   Display data envelopment analysis results stored in a 'deaout'
#   structure.
#   DEADISP( out ) Display data envelopment analysis results.
#   DEADISP( out, dispstr ) Display results using the specified 'dispstr'.
    
    #   Example
#       
#       io = dea(X, Y, 'orient', 'io');
#       deadisp(io);
    
    #       deadisp(io, 'names/lambda/eff');
    
    #   See also DEAOUT, DEA2TABLE
    
    #   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
#   http://www.deatoolbox.com
    
    #   Version: 1.0
#   LAST UPDATE: 9, May, 2017
    
    # Check if input is a structure
    if logical_not(isstruct(out)):
        error('Input is not a structure')
    
    
    # If not custom dispstr is specified
    if nargin < 2:
        # Check if the structure has a dispstr
        if isfield(out,'dispstr'):
            dispstr=out.dispstr
# .\deadisp.m:33
        else:
            error('Input structure does not have a display string, \'dispstr\', field')
    
    # TITLE
    fprintf('_______________________________\\n')
    if isfield(out,'disptext_title'):
        # Display specified title
        fprintf('<strong>%s</strong>\\n\\n',out.disptext_title)
    else:
        fprintf('<strong>Data Envelopment Analysis (DEA)</strong>\\n\\n')
    
    
    # TEXT 1: Before model information
    if isfield(out,'disptext_text1'):
        fprintf('%s\\n',out.disptext_text1)
    
    
    # MODEL INFORMATION    
    # DMU and number of inputs and outputs information
    if isfield(out,'n'):
        fprintf('DMUs: %i ',out.n)
        if isfield(out,'neval'):
            if out.neval != out.n and logical_not(isnan(out.neval)):
                fprintf('(%i evaluated)',out.neval)
        fprintf('\\n')
    
    
    if isfield(out,'m') and isfield(out,'s'):
        fprintf('Inputs: %i     Outputs: %i ',out.m,out.s)
        if isfield(out,'r'):
            if logical_not(isnan(out.r)):
                fprintf('    Undesirable: %i ',out.r)
        fprintf('\\n')
    
    
    # Model
    if isfield(out,'model'):
        fprintf('Model: %s ',out.model)
        fprintf('\\n')
    
    
    # Orientation
    if isfield(out,'orient'):
        fprintf('Orientation: %s ',out.orient)
        if 'io' == (out.orient):
            fprintf('(Input oriented)')
        else:
            if 'oo' == (out.orient):
                fprintf('(Output oriented)')
            else:
                if 'ddf' == (out.orient):
                    fprintf('(Directional distance function)')
        fprintf('\\n')
    
    
    # Returns to scale
    if isfield(out,'rts'):
        fprintf('Returns to scale: %s ',out.rts)
        if 'crs' == (out.rts):
            fprintf('(Constant)')
        else:
            if 'vrs' == (out.rts):
                fprintf('(Variable)')
            else:
                if 'scaleeff' == (out.rts):
                    fprintf('(Scale efficiency)')
        fprintf('\\n')
    
    
    # Bootstrap and significance
    if isfield(out,'nreps'):
        if logical_not(isnan(out.nreps)):
            fprintf('Bootstrap replications: %i \\n',out.nreps)
    
    if isfield(out,'alpha'):
        if logical_not(isnan(out.alpha)):
            fprintf('Significance level: %4.2f \\n',out.alpha)
    
    
    fprintf('\\n')
    
    if isfield(out,'disptext_text2'):
        fprintf('%s\\n',out.disptext_text2)
    
    
    # Period (for temporal models)
    if isfield(out,'period'):
        if 'base' == (out.period):
            disp('Reference period is base period')
        else:
            if 'comparison' == (out.period):
                disp('Reference period is comparison period')
            else:
                if 'geomean' == (out.period):
                    disp('Geometric mean is computed')
    
    
    # Fixbase t (for temporal models)
    if isfield(out,'fixbaset'):
        if out.fixbaset == 1:
            disp('Base period is period 1')
        else:
            disp('Base period is previous period')
        disp(' ')
    
    
    # TEXT 3: Before table
    if isfield(out,'disptext_text3'):
        fprintf('%s\\n',out.disptext_text3)
    
    
    # TABLE
    dispstr=strsplit(dispstr,'/')
# .\deadisp.m:156
    tabAll=[]
# .\deadisp.m:157
    for i in arange(1,length(dispstr)).reshape(-1):
        # Get param name
        paramstr=char(dispstr(i))
# .\deadisp.m:160
        name,format=getDEAformat(paramstr,out.orient,nargout=2)
# .\deadisp.m:163
        if isempty(name):
            disptext_field=sprintf('disptext_%s',strrep(paramstr,'.','_'))
# .\deadisp.m:166
            if isfield(out,disptext_field):
                # If custom name exists in the output structure use it
                name=eval(sprintf('out.%s',disptext_field))
# .\deadisp.m:169
            else:
                # If not, display paramstr name without eff.
                name=strrep(paramstr,'eff.','')
# .\deadisp.m:172
        # Get data
        dat=eval(sprintf('out.%s',paramstr))
# .\deadisp.m:177
        if logical_not(iscell(dat)):
            # Convert to cell if not cell
            dat=num2cell(dat)
# .\deadisp.m:180
        # Number of columns
        ncols=size(dat,2)
# .\deadisp.m:184
        for j in arange(1,ncols).reshape(-1):
            # Get Body
            bodyj=cellfun(lambda x=None: sprintf(format,x),dat(arange(),j),'Unif',false)
# .\deadisp.m:190
            if ncols > 1:
                # If more than 1 columns add number to name
                namej=concat([name,num2str(j)])
# .\deadisp.m:195
            else:
                namej=copy(name)
# .\deadisp.m:197
            headerj=cellstr(namej)
# .\deadisp.m:200
            allThis_c=concat([[headerj],[' '],[bodyj]])
# .\deadisp.m:203
            tabThis=char(allThis_c)
# .\deadisp.m:206
            tabThis=strjust(tabThis,'right')
# .\deadisp.m:209
            tabThis[arange(),arange(2,end() + 1)]=tabThis(arange(),arange())
# .\deadisp.m:212
            tabThis[arange(),1]=' '
# .\deadisp.m:213
            tabThis[arange(),end() + 1]='|'
# .\deadisp.m:214
            # Add to table
            tabAll=concat([tabAll,tabThis])
# .\deadisp.m:217
        # Replace second row with separator
        tabAll[2,arange()]='-'
# .\deadisp.m:222
    
    
    disp(repelem('-',size(tabAll,2)))
    
    disp(tabAll)
    
    disp(repelem('-',size(tabAll,2)))
    
    # TEXT 4: After table
    if isfield(out,'disptext_text4'):
        fprintf('%s\\n',out.disptext_text4)
    
    
    
    return
    
if __name__ == '__main__':
    pass
    