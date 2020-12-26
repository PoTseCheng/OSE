import pandas as pd

def reading(x):
    '''
    reading the csv files
    '''
    df=pd.read_csv(x)
    add=df.columns.tolist()
    df.loc[-1] = add 
    df.index = df.index + 1  
    df = df.sort_index()
    df.columns=[x.replace(".csv","")]

    return df


