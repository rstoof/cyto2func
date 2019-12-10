import csv
import numpy as np
import matplotlib.pyplot as plt
import re
import os
import multiprocessing
import pandas
import fcsparser
import seaborn
import pandas
from datetime import datetime
import scipy.stats as st
%matplotlib inline


import numpy as np

####
df=pandas.read_csv('file_description.csv')
df.filename=df.filename.replace({'.mqd':'.fcs'}, regex=True)
df.filename2=df.filename.replace({'.mqd':'.fcs'}, regex=True)

channels=["GFP_H"]

median_fcs=[]
droparr=[]
for index,row in df.head().iterrows():
    try:
        meta, data = fcsparser.parse("/home/ruud/Desktop/server/home/data/Local_Limited_storage/FCS/"+row.filename, meta_data_only=False, reformat_meta=True)
        data.columns=[x.strip().replace('-', '_') for x in data.columns]
        medi=data["GFP_H"].median()
        data["GFP_H"].hist()
        print(medi)
        median_fcs.append(medi)
    except:
        print("fail at"+str(index))
        droparr.append(index)
dfjm=df.copy()
dfjm=dfjm.drop(droparr)
dfjm
dfjm.insert(6,"median",median_fcs)
dfjm.to_csv("just_medians.csv")
