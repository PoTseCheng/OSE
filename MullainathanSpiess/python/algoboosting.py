import xgboost as xgb
import time

parameter = list(eta = .3, max_depth = 6, nround = 1000, gamma=0)
def xgbcontinousalg(df, LHS, RHS, f, to, quiet = True):
    ''''''

    if quiet != True:
        start = time.time()
    parameter={ "eta" : [.3],
                "max_depth" : [6], 
                "nround" : [1000], 
                "gamma":[0]}
    formula=[str(LHS)+"~ "+ list(set(RHS)-set(parameters["interactions"]))]

