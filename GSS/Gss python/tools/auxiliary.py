import pandas as pd
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def reading(x):
    '''
    Homemade reader for csv file produced by homemade converter, an alternative for scipy.io.loadmat.

    Parameters:
    x(str):Csv file from converter

    Returns:
    Pandas dataframe
    '''
    df=pd.read_csv(x)
    add=df.columns.tolist()
    df.loc[-1] = add 
    df.index = df.index + 1  
    df = df.sort_index()
    df.columns=[x.replace(".csv","")]

    return df

####Figures########################################################

def Figure1():
    '''
    This function showcase the ergodic sets
    '''
    df1 = pd.read_csv("country_k.csv")
    df2 = pd.read_csv("aT20200N10.csv")
    df1.drop(df1.tail(1).index,inplace=True)
    df2=df2.loc[ :,:"country_2"]
    df2= df2.head(2000)
    result1 = pd.concat([df1["country_1"], df2["country_1"]], axis=1)
    result1.columns=["Capital", "Productivity"]
    result1["Country"]="country 1"
    Fig1= sns.jointplot("Capital", "Productivity",data=result1, kind='hex')

    return
