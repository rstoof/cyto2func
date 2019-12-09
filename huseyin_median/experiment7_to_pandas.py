import pandas
import os
from tkinter import *
from tkinter import filedialog


root = Tk()
root.withdraw()
yourdir = filedialog.askdirectory(title = "Select Measurement folder")
usedir=yourdir

file_list=os.listdir(usedir)
file_list
df_list=[]
for file in file_list:
    df=pandas.read_csv(usedir+"/"+file,header=None,names=["iptg","median_yfp"])
    df["filename"]=file
    df_list.append(df)
big_df = pandas.concat(df_list)

big_df["strain"]=big_df.filename.str.extract(r"op_(.*?)_")
big_df["backbone"]=big_df.filename.str.extract(r"_.*?_(.*?)_")
big_df["plasmid"]=big_df.filename.str.extract(r"_.*?_.*?_(.*)?\.")
big_df
big_df.to_csv("huseyin_median.csv")
