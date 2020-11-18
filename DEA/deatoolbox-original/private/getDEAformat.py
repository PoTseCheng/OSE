# Generated with SMOP  0.41
from libsmop import *
# .\getDEAformat.m

    
@function
def getDEAformat(fieldname=None,orient=None,*args,**kwargs):
    varargin = getDEAformat.varargin
    nargin = getDEAformat.nargin

    #GETDEAFORMAT Private function
#   Private function
    
    #   Copyright 2016 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
#   http://www.deatoolbox.com
    
    #   Version: 1.0
#   LAST UPDATE: 6, May, 2017
    
    if nargin < 2:
        orient='none'
# .\getDEAformat.m:13
    
    
    fmtNumber='%7.4f'
# .\getDEAformat.m:16
    if 'names' == fieldname:
        name='DMU'
# .\getDEAformat.m:21
        format='%s'
# .\getDEAformat.m:22
    else:
        if 'X' == fieldname:
            name='X'
# .\getDEAformat.m:24
            format=copy(fmtNumber)
# .\getDEAformat.m:25
        else:
            if 'Y' == fieldname:
                name='Y'
# .\getDEAformat.m:27
                format=copy(fmtNumber)
# .\getDEAformat.m:28
            else:
                if 'eff' == fieldname:
                    if 'io' == orient:
                        name='Theta'
# .\getDEAformat.m:32
                    else:
                        if 'oo' == orient:
                            name='Phi'
# .\getDEAformat.m:34
                        else:
                            if 'ddf' == orient:
                                name='Beta'
# .\getDEAformat.m:36
                            else:
                                name='Eff'
# .\getDEAformat.m:38
                    format=fmtNumber.T
# .\getDEAformat.m:40
                else:
                    if 'slack.X' == fieldname:
                        name='slackX'
# .\getDEAformat.m:42
                        format=copy(fmtNumber)
# .\getDEAformat.m:43
                    else:
                        if 'slack.Y' == fieldname:
                            name='slackY'
# .\getDEAformat.m:45
                            format=copy(fmtNumber)
# .\getDEAformat.m:46
                        else:
                            if 'lambda' == fieldname:
                                name='lambda'
# .\getDEAformat.m:48
                                format=copy(fmtNumber)
# .\getDEAformat.m:49
                            else:
                                if 'Xeff' == fieldname:
                                    name='Xeff'
# .\getDEAformat.m:51
                                    format=copy(fmtNumber)
# .\getDEAformat.m:52
                                else:
                                    if 'Yeff' == fieldname:
                                        name='Yeff'
# .\getDEAformat.m:54
                                        format=copy(fmtNumber)
# .\getDEAformat.m:55
                                    else:
                                        if 'dual.X' == fieldname:
                                            name='dualX'
# .\getDEAformat.m:57
                                            format=copy(fmtNumber)
# .\getDEAformat.m:58
                                        else:
                                            if 'dual.Y' == fieldname:
                                                name='dualY'
# .\getDEAformat.m:60
                                                format=copy(fmtNumber)
# .\getDEAformat.m:61
                                            else:
                                                if 'dual.rts' == fieldname:
                                                    name='dualRTS'
# .\getDEAformat.m:63
                                                    format=copy(fmtNumber)
# .\getDEAformat.m:64
                                                else:
                                                    if 'exitflag' == fieldname:
                                                        name='EFlag'
# .\getDEAformat.m:66
                                                        format='%i'
# .\getDEAformat.m:67
                                                    else:
                                                        if 'eff.crs' == fieldname:
                                                            name='CRS'
# .\getDEAformat.m:70
                                                            format=copy(fmtNumber)
# .\getDEAformat.m:71
                                                        else:
                                                            if 'eff.vrs' == fieldname:
                                                                name='VRS'
# .\getDEAformat.m:73
                                                                format=copy(fmtNumber)
# .\getDEAformat.m:74
                                                            else:
                                                                if 'eff.scale' == fieldname:
                                                                    name='ScaleEff'
# .\getDEAformat.m:76
                                                                    format=copy(fmtNumber)
# .\getDEAformat.m:77
                                                                else:
                                                                    if 'eff.M' == fieldname:
                                                                        name='M'
# .\getDEAformat.m:80
                                                                        format=copy(fmtNumber)
# .\getDEAformat.m:81
                                                                    else:
                                                                        if 'eff.MTEC' == fieldname:
                                                                            name='MTEC'
# .\getDEAformat.m:83
                                                                            format=copy(fmtNumber)
# .\getDEAformat.m:84
                                                                        else:
                                                                            if 'eff.MTC' == fieldname:
                                                                                name='MTC'
# .\getDEAformat.m:86
                                                                                format=copy(fmtNumber)
# .\getDEAformat.m:87
                                                                            else:
                                                                                if 'Xprice' == fieldname:
                                                                                    name='Xprice'
# .\getDEAformat.m:90
                                                                                    format=copy(fmtNumber)
# .\getDEAformat.m:91
                                                                                else:
                                                                                    if 'Yprice' == fieldname:
                                                                                        name='Yprice'
# .\getDEAformat.m:93
                                                                                        format=copy(fmtNumber)
# .\getDEAformat.m:94
                                                                                    else:
                                                                                        if 'eff.C' == fieldname:
                                                                                            name='CostEff'
# .\getDEAformat.m:96
                                                                                            format=copy(fmtNumber)
# .\getDEAformat.m:97
                                                                                        else:
                                                                                            if 'eff.R' == fieldname:
                                                                                                name='RevEff'
# .\getDEAformat.m:99
                                                                                                format=copy(fmtNumber)
# .\getDEAformat.m:100
                                                                                            else:
                                                                                                if 'eff.P' == fieldname:
                                                                                                    name='ProfEff'
# .\getDEAformat.m:102
                                                                                                    format=copy(fmtNumber)
# .\getDEAformat.m:103
                                                                                                else:
                                                                                                    if 'eff.A' == fieldname:
                                                                                                        name='AllocEff'
# .\getDEAformat.m:105
                                                                                                        format=copy(fmtNumber)
# .\getDEAformat.m:106
                                                                                                    else:
                                                                                                        if 'eff.T' == fieldname:
                                                                                                            name='TechEff'
# .\getDEAformat.m:108
                                                                                                            format=copy(fmtNumber)
# .\getDEAformat.m:109
                                                                                                        else:
                                                                                                            if 'Yu' == fieldname:
                                                                                                                name='Yu'
# .\getDEAformat.m:112
                                                                                                                format=copy(fmtNumber)
# .\getDEAformat.m:113
                                                                                                            else:
                                                                                                                if 'slack.Yu' == fieldname:
                                                                                                                    name='slackYu'
# .\getDEAformat.m:115
                                                                                                                    format=copy(fmtNumber)
# .\getDEAformat.m:116
                                                                                                                else:
                                                                                                                    if 'Yueff' == fieldname:
                                                                                                                        name='Yueff'
# .\getDEAformat.m:118
                                                                                                                        format=copy(fmtNumber)
# .\getDEAformat.m:119
                                                                                                                    else:
                                                                                                                        if 'eff.ML' == fieldname:
                                                                                                                            name='ML'
# .\getDEAformat.m:122
                                                                                                                            format=copy(fmtNumber)
# .\getDEAformat.m:123
                                                                                                                        else:
                                                                                                                            if 'eff.MLTEC' == fieldname:
                                                                                                                                name='MLTEC'
# .\getDEAformat.m:125
                                                                                                                                format=copy(fmtNumber)
# .\getDEAformat.m:126
                                                                                                                            else:
                                                                                                                                if 'eff.MLTC' == fieldname:
                                                                                                                                    name='MLTC'
# .\getDEAformat.m:128
                                                                                                                                    format=copy(fmtNumber)
# .\getDEAformat.m:129
                                                                                                                                else:
                                                                                                                                    if 'eff.o' == fieldname:
                                                                                                                                        name='eff'
# .\getDEAformat.m:132
                                                                                                                                        format=copy(fmtNumber)
# .\getDEAformat.m:133
                                                                                                                                    else:
                                                                                                                                        if 'eff.b' == fieldname:
                                                                                                                                            name='effboot'
# .\getDEAformat.m:135
                                                                                                                                            format=copy(fmtNumber)
# .\getDEAformat.m:136
                                                                                                                                        else:
                                                                                                                                            if 'eff.c' == fieldname:
                                                                                                                                                name='effCI'
# .\getDEAformat.m:138
                                                                                                                                                format=copy(fmtNumber)
# .\getDEAformat.m:139
                                                                                                                                            else:
                                                                                                                                                if 'eff.M.o' == fieldname:
                                                                                                                                                    name='M'
# .\getDEAformat.m:142
                                                                                                                                                    format=copy(fmtNumber)
# .\getDEAformat.m:143
                                                                                                                                                else:
                                                                                                                                                    if 'eff.M.b' == fieldname:
                                                                                                                                                        name='Mboot'
# .\getDEAformat.m:145
                                                                                                                                                        format=copy(fmtNumber)
# .\getDEAformat.m:146
                                                                                                                                                    else:
                                                                                                                                                        if 'eff.M.cL' == fieldname:
                                                                                                                                                            name='McLow'
# .\getDEAformat.m:148
                                                                                                                                                            format=copy(fmtNumber)
# .\getDEAformat.m:149
                                                                                                                                                        else:
                                                                                                                                                            if 'eff.M.cU' == fieldname:
                                                                                                                                                                name='McUpp'
# .\getDEAformat.m:151
                                                                                                                                                                format=copy(fmtNumber)
# .\getDEAformat.m:152
                                                                                                                                                            else:
                                                                                                                                                                if 'eff.MTEC.o' == fieldname:
                                                                                                                                                                    name='MTEC'
# .\getDEAformat.m:154
                                                                                                                                                                    format=copy(fmtNumber)
# .\getDEAformat.m:155
                                                                                                                                                                else:
                                                                                                                                                                    if 'eff.MTEC.b' == fieldname:
                                                                                                                                                                        name='MTECboot'
# .\getDEAformat.m:157
                                                                                                                                                                        format=copy(fmtNumber)
# .\getDEAformat.m:158
                                                                                                                                                                    else:
                                                                                                                                                                        if 'eff.MTEC.cL' == fieldname:
                                                                                                                                                                            name='MTECcLow'
# .\getDEAformat.m:160
                                                                                                                                                                            format=copy(fmtNumber)
# .\getDEAformat.m:161
                                                                                                                                                                        else:
                                                                                                                                                                            if 'eff.MTEC.cU' == fieldname:
                                                                                                                                                                                name='MTECcUpp'
# .\getDEAformat.m:163
                                                                                                                                                                                format=copy(fmtNumber)
# .\getDEAformat.m:164
                                                                                                                                                                            else:
                                                                                                                                                                                if 'eff.MTC.o' == fieldname:
                                                                                                                                                                                    name='MTC'
# .\getDEAformat.m:166
                                                                                                                                                                                    format=copy(fmtNumber)
# .\getDEAformat.m:167
                                                                                                                                                                                else:
                                                                                                                                                                                    if 'eff.MTC.b' == fieldname:
                                                                                                                                                                                        name='MTCboot'
# .\getDEAformat.m:169
                                                                                                                                                                                        format=copy(fmtNumber)
# .\getDEAformat.m:170
                                                                                                                                                                                    else:
                                                                                                                                                                                        if 'eff.MTC.cL' == fieldname:
                                                                                                                                                                                            name='MTCcLow'
# .\getDEAformat.m:172
                                                                                                                                                                                            format=copy(fmtNumber)
# .\getDEAformat.m:173
                                                                                                                                                                                        else:
                                                                                                                                                                                            if 'eff.MTC.cU' == fieldname:
                                                                                                                                                                                                name='MTCcUpp'
# .\getDEAformat.m:175
                                                                                                                                                                                                format=copy(fmtNumber)
# .\getDEAformat.m:176
                                                                                                                                                                                            else:
                                                                                                                                                                                                name=[]
# .\getDEAformat.m:178
                                                                                                                                                                                                format=copy(fmtNumber)
# .\getDEAformat.m:179
    
    
    return name,format
    
if __name__ == '__main__':
    pass
    