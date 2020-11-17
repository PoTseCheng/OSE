# Generated with SMOP  0.41
from libsmop import *
# .\deacross.m

    
@function
def deacross(X=None,Y=None,varargin=None,*args,**kwargs):
    varargin = deacross.varargin
    nargin = deacross.nargin

    #DEACROSS Data envelopment analysis cross efficiency model
#   Computes data envelopment analysis cross efficiency model using
#   Sexton's et al. (1986) model.
    
    #   out = DEACROSS(X, Y, Name, Value) computes data envelopment analysis 
#   cross efficiency model with inputs X and outputs Y. Model properties are
#   specified using one or more Name ,Value pair arguments.
    
    #   Additional properties:
#   - 'orient': orientation. Input oriented 'io' (Default).
#   - 'model': 'linear' (Default).
#   - 'objective': minimization of the objective function, aggresive
#   approach ('aggresive'). Maximization of the objective function,
#   benevolent approach ('benevolent').
#   - 'mean': include the evaluated DMU when computing the cross efficiency 
#   score ('inclusive'); or excludes it ('exclusive'). Geometric mean 
#   including the evaluated DMU ('ginclusive'); or geometric mean excluding
#   it ('gexclusive').
#   - 'names': DMU names.
    
    
    #   Example
#     
#      cross_agg = deacross(X, Y, 'objective', 'aggressive', 'mean', 'inclusive');
#      deadisp(cross_agg);
    
    #      cross_ben = deacross(X, Y, 'objective', 'benevolent', 'mean', 'exclusive');
#      deadisp(cross_ben);
    
    #   See also DEAOUT, DEA, DEABOOT
    
    #   Copyright 2018 Inmaculada C. Alvarez, Javier Barbero, Jose L. Zofio
#   http://www.deatoolbox.com
    
    #   Version: 1.0
#   LAST UPDATE: 10, July, 2018
    
    # Check size
    if size(X,1) != size(Y,1):
        error('Number of rows in X must be equal to number of rows in Y')
    
    
    # Get number of DMUs (n), inputs (m) and outputs (s)
    n,m=size(X,nargout=2)
# .\deacross.m:46
    s=size(Y,2)
# .\deacross.m:47
    
    # Default optimization-options
    optimoptsdef=optimoptions('linprog','display','off','Algorithm','dual-simplex','TolFun',1e-10,'TolCon',1e-07)
# .\deacross.m:51
    
    p=copy(inputParser)
# .\deacross.m:54
    addPar=lambda v1=None,v2=None,v3=None,v4=None: addParameter(v1,v2,v3,v4)
# .\deacross.m:55
    
    addPar(p,'names',cellstr(int2str((arange(1,n)).T)),lambda x=None: iscellstr(x) and (length(x) == n))
    addPar(p,'optimopts',optimoptsdef,lambda x=None: logical_not(isempty(x)))
    
    addPar(p,'orient','io',lambda x=None: any(validatestring(x,cellarray(['io','oo']))))
    addPar(p,'model','linear',lambda x=None: any(validatestring(x,cellarray(['linear']))))
    addPar(p,'objective','aggressive',lambda x=None: any(validatestring(x,cellarray(['aggressive','benevolent']))))
    addPar(p,'mean','inclusive',lambda x=None: any(validatestring(x,cellarray(['inclusive','exclusive','ginclusive','gexclusive']))))
    p.parse(varargin[arange()])
    options=p.Results
# .\deacross.m:71
    
    if size(options.names,2) > 1:
        options.names = copy(options.names.T)
# .\deacross.m:75
    
    
    
    if strcmp(options.orient,'oo'):
        error('Output oriented orientation not yet implemented')
    
    
    # OPTIMIZATION OPTIONS:
    optimopts=options.optimopts
# .\deacross.m:84
    
    if 'aggressive' == (options.objective):
        appSign=1
# .\deacross.m:89
    else:
        if 'benevolent' == (options.objective):
            appSign=- 1
# .\deacross.m:91
    
    # LINEAR APPROACH
    if 'linear' == (options.model):
        # INPUT-OTINETED DEA
        eff=nan(n,1)
# .\deacross.m:98
        v=nan(n,m)
# .\deacross.m:99
        u=nan(n,s)
# .\deacross.m:100
        for j in arange(1,n).reshape(-1):
            # Objective function
            f=- concat([zeros(1,m),Y(j,arange())])
# .\deacross.m:106
            A=concat([- X,Y])
# .\deacross.m:109
            b=concat([zeros(n,1)])
# .\deacross.m:110
            Aeq=concat([X(j,arange()),zeros(1,s)])
# .\deacross.m:111
            beq=1
# .\deacross.m:112
            lb=zeros(1,m + s)
# .\deacross.m:113
            z,__,exitflag=linprog(f,A,b,Aeq,beq,lb,[],[],optimopts,nargout=3)
# .\deacross.m:116
            if exitflag != 1:
                if options.warning:
                    warning('DMU %i. First Step. Optimization exit flag: %i',j,exitflag)
            # Store efficiency and weights
            v[j,arange()]=z(arange(1,m))
# .\deacross.m:124
            u[j,arange()]=z(arange(m + 1,m + s))
# .\deacross.m:125
            eff[j]=sum(multiply(z(arange(m + 1,m + s)),Y(j,arange()).T))
# .\deacross.m:126
        # AGGREGATE INPUTS AND OUTPUTs
        Xagg=sum(X) - X
# .\deacross.m:132
        Yagg=sum(Y) - Y
# .\deacross.m:133
        v_ab=nan(n,m)
# .\deacross.m:136
        u_ab=nan(n,s)
# .\deacross.m:137
        for j in arange(1,n).reshape(-1):
            # Objective function
            f=dot(appSign,concat([- Xagg(j,arange()),Yagg(j,arange())]))
# .\deacross.m:142
            A=concat([- X,Y])
# .\deacross.m:145
            b=concat([zeros(n,1)])
# .\deacross.m:146
            Aeq=concat([[X(j,arange()),zeros(1,s)],[dot(- eff(j),X(j,arange())),Y(j,arange())]])
# .\deacross.m:147
            beq=concat([[1],[0]])
# .\deacross.m:149
            lb=zeros(1,m + s)
# .\deacross.m:151
            z,__,exitflag=linprog(f,A,b,Aeq,beq,lb,[],[],optimopts,nargout=3)
# .\deacross.m:154
            if exitflag != 1:
                if options.warning:
                    warning('DMU %i. First Step. Optimization exit flag: %i',j,exitflag)
            # Store weights
            v_ab[j,arange()]=z(arange(1,m))
# .\deacross.m:162
            u_ab[j,arange()]=z(arange(m + 1,m + s))
# .\deacross.m:163
        # Peer-Appraisal
        PA=(dot(u_ab,Y.T)) / (dot(v_ab,X.T))
# .\deacross.m:168
        PAex=copy(PA)
# .\deacross.m:171
        PAex=reshape(PAex,dot(n,n),1)
# .\deacross.m:172
        PAex[arange(1,dot(n,n),n + 1)]=[]
# .\deacross.m:173
        PAex=reshape(PAex,n - 1,n)
# .\deacross.m:174
        if 'inclusive' == (options.mean):
            crosseff=mean(PA).T
# .\deacross.m:178
        else:
            if 'exclusive' == (options.mean):
                crosseff=mean(PAex).T
# .\deacross.m:180
            else:
                if 'ginclusive' == (options.mean):
                    crosseff=geomean(PA).T
# .\deacross.m:182
                else:
                    if 'gexclusive' == (options.mean):
                        crosseff=geomean(PAex).T
# .\deacross.m:184
    
    
    # Generate output structure
    out.n = copy(n)
# .\deacross.m:190
    out.m = copy(m)
# .\deacross.m:191
    out.s = copy(s)
# .\deacross.m:192
    out.names = copy(options.names)
# .\deacross.m:194
    out.model = copy('deacross-linear')
# .\deacross.m:196
    out.orient = copy(options.orient)
# .\deacross.m:197
    out.rts = copy('crs')
# .\deacross.m:198
    
    effs.eff = copy(eff)
# .\deacross.m:201
    effs.crosseff = copy(crosseff)
# .\deacross.m:202
    effs.PA = copy(PA)
# .\deacross.m:203
    effs.vab = copy(v_ab)
# .\deacross.m:204
    effs.uab = copy(u_ab)
# .\deacross.m:205
    
    out.eff = copy(effs)
# .\deacross.m:208
    
    out.dispstr = copy('names/eff.eff/eff.crosseff')
# .\deacross.m:211
    
    out.disptext_title = copy('Data Envelopment Analysis (DEA)')
# .\deacross.m:214
    out.disptext_text2 = copy('Cross efficiency')
# .\deacross.m:215
    if 'aggressive' == (options.objective):
        out.disptext_text2 = copy(concat([out.disptext_text2,': Aggressive Approach']))
# .\deacross.m:218
    else:
        if 'benevolent' == (options.objective):
            out.disptext_text2 = copy(concat([out.disptext_text2,': Benevolent Approach']))
# .\deacross.m:220
    
    
    out.disptext_eff_eff = copy('Eff')
# .\deacross.m:223
    out.disptext_eff_crosseff = copy('CrossEff')
# .\deacross.m:224
    out.disptext_eff_PA = copy('PA')
# .\deacross.m:225
    out.disptext_eff_vab = copy('v_ab')
# .\deacross.m:226
    out.disptext_eff_uab = copy('u_ab')
# .\deacross.m:227
    out.disptext_text4 = copy('CrossEff = Cross Efficiency')
# .\deacross.m:229
    if 'inclusive' == (options.mean):
        out.disptext_text4 = copy(concat([out.disptext_text4,' (inclusive mean)']))
# .\deacross.m:232
    else:
        if 'exclusive' == (options.mean):
            out.disptext_text4 = copy(concat([out.disptext_text4,' (exclusive mean)']))
# .\deacross.m:234
        else:
            if 'ginclusive' == (options.mean):
                out.disptext_text4 = copy(concat([out.disptext_text4,' (inclusive geometric mean)']))
# .\deacross.m:236
            else:
                if 'gexclusive' == (options.mean):
                    out.disptext_text4 = copy(concat([out.disptext_text4,' (exclusive geometric mean)']))
# .\deacross.m:238
    
    
    
    return out
    
if __name__ == '__main__':
    pass
    