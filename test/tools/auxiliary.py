# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 15:38:26 2020

@author: Viktor Cheng
"""
#I did not coded everything by myself!

import os


#####

def get_data_storage():
    """
    :return: The absolute path of the folder data.
    """
    dirname = os.path.dirname(__file__)
    path = dirname + "\pkl\group_data"
    return path



