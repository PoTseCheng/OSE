

####this will be the python equivalent for dea.m



def dea(X, Y, varargin):
    '''
    python equivalent for dea.m
    '''
    if X.shape[0] != Y.shape[0]:
        error='Number of rows in X must be equal to number of rows in Y'
        return error
    
    
    
    options = getDEAoptions(n, varargin{:})