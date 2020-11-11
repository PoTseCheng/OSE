# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:15:12 2020

@author: Viktor Cheng
"""

#This python file is a simplified version of the original file from Thomas H. Joegensen,
#the purpose is to combine the existing code with the ruspy package developed by the OSE initiative.

import pandas as pd
from numba import njit, jitclass, prange, boolean, int32, double #blah blah
import numpy as np
import matplotlib.pyplot as plot
import ruspy as rpy
import seaborn as sns



#####Preparations#####

#One should also use the csv file that I constructed out of the original model

def read(x):
    """
    The read function reads in the csv file and automatically unpack it 
    
    
    """
    df=pd.read_csv("x")
    mom_data= np.log(df.iloc[:,0].to_numpy())
    weight= df.iloc[:,2:43].to_numpy()
    income= df.iloc[:1].to_numpy()
    theta0 = [0.944, 1.860] # need to keep this line
    
    return mom_data, weight, income, theta0
    



