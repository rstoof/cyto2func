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

df=pandas.read_csv('file_description.csv')
df.filename=df.filename.replace({'.mqd':'.fcs'}, regex=True)
df.filename2=df.filename.replace({'.mqd':'.fcs'}, regex=True)

channels=["FSC_H","GFP_H"]

#######functions
def kdedata(dat,usestandardrange=False): #makes a 2D histogram of log transformed data
    [fl,vl]=dat
    if usestandardrange:
        fluormin =np.log(.001)
        fluormax =np.log(100)
    else:
        fluormin =np.percentile(fl,5)
        fluormax = np.percentile(fl,95)
    volmin = np.percentile(vl,5)
    volmax =np.percentile(vl,95)
    #dat=dat[dat[:,0]>fluormin&dat[:,0]<fluormax&dat[:,1]>volmin&dat[:,1]<volmax]
    xx, yy = np.mgrid[fluormin:fluormax:160j, volmin:volmax:160j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([fl, vl])
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    #plt.imshow(np.rot90(f),aspect="auto", cmap=plt.cm.gist_earth_r,extent=[fluormin, fluormax, volmin, volmax])
    #plt.show()
    return [xx,yy,f]
def fit2binormal(kdedat):
    #fits a 2D -normal distribution to 2D histogram of log transformed data
    [xx,yy,f]=kdedat
    xy_mesh = [xx,yy]
    guess_vals = [np.max(f), 3, 3, 2, 2, .8]
    #guess_vals = [np.max(f), 2, 3, 1, 1, .6]
    # perform the fit, making sure to flatten the noisy data for the fit routine
    fit_params, cov_mat = curve_fit(gaussian_2d, xy_mesh, np.ravel(f), p0=guess_vals, maxfev=10000000,bounds=([0,-3,-3,0,0,.0],[50,15,15,5,5,1]))
    # calculate fit parameter errors from covariance matrix
    fit_errors = np.sqrt(np.diag(cov_mat))
    # manually calculate R-squared goodness of fit
    fit_residual = f - gaussian_2d(xy_mesh, *fit_params).reshape(np.outer(xx[:,0],yy[0]).shape)
    fit_Rsquared = 1 - np.var(fit_residual)/np.var(f)
    return np.append(fit_params,fit_Rsquared)
from scipy.optimize import curve_fit
def gaussian_2d(xy_mesh, amp, xc, yc, sigma_x, sigma_y,rho):
    # unpack 1D list into 2D x and y coords
    (x, y) = xy_mesh
    # make the 2D Gaussian matrix
    gauss = amp*np.exp(-((x-xc)**2/(sigma_x**2)-(2*rho*(x-xc)*(y-yc))/(sigma_x*sigma_y)+(y-yc)**2/(sigma_y**2)))/(2*np.pi*sigma_x*sigma_y*np.sqrt(1-rho**2))
    # flatten the 2D Gaussian down to 1D
    return np.ravel(gauss)
########code
df=pandas.read_csv('file_description.csv')
df.filename=df.filename.replace({'.mqd':'.fcs'}, regex=True)
#df.filename2=df.filename.replace({'.mqd':'.fcs'}, regex=True)

channels=["FSC_H","GFP_H"]
testfile=df.head(1).filename[0]
testmeta,testdata=fcsparser.parse("../FCS/"+testfile, meta_data_only=False, reformat_meta=True)
testdata.columns=[x.strip().replace('-', '_') for x in testdata.columns]
testdata=testdata.query(channels[0]+">0 and "+channels[1]+">0")
testdata2=kdedata([np.log(testdata["GFP_H"]),np.log(testdata["FSC_H"])])
fit2binormal(testdata2)

minvalfsc=[]
minvalgfp=[]
datearr=[]
droparr=[]
fitarr=[]
gate=True

for index,row in df.iterrows():
    try:
        meta, data = fcsparser.parse("../FCS/"+row.filename, meta_data_only=False, reformat_meta=True)
        data.columns=[x.strip().replace('-', '_') for x in data.columns]
        if gate==True:
            data=data[(data["SSC_H"]>np.exp(2.5)) & (data["SSC_A"]>0)&(data["FSC_H"]>np.exp(1.5))&(data["GFP_H"]>0)&(data["SSC_A"]>data["SSC_H"])&(data["SSC_H"]>data["FSC_H"])]
        datetime_object = datetime.strptime(meta['$DATE']+" "+meta['$BTIM'], '%Y-%b-%d %H:%M:%S')
        datearr.append(datetime_object)
        minvalfsc.append(data.SSC_H.quantile(.05))
        minvalgfp.append(data.GFP_H.quantile(.05))
        data=data.query(channels[0]+">0 and "+channels[1]+">0")
        data2=kdedata([np.log(data["GFP_H"]),np.log(data["FSC_H"])])
        fits=fit2binormal(data2)
        fitarr.append(fits)
        print(index)
        print(fits)
    except:
        print("fail at"+str(index))
        droparr.append(index)
df
minvalfsc
type(fitarr)
df3=df.copy()
df3=df3.drop(droparr)
df3
df3.insert(6,"lowfsc",minvalfsc)
df3.insert(7,"lowgfp",minvalgfp)
df3.insert(8,"real_time",datearr)

fitmufl=[fit[1] for fit in fitarr]
fitmuv=[fit[2] for fit in fitarr]
fitstdfl=[fit[3] for fit in fitarr]
fitstdv=[fit[4] for fit in fitarr]
fitrho=[fit[5] for fit in fitarr]
fitgoodness=[fit[6] for fit in fitarr]
fitmufl[1]

os.getcwd()
df3.insert(9,"log_mean_gfp",fitmufl)
df3.insert(10,"log_mean_v",fitmuv)
df3.insert(11,"log_std_gfp",fitstdfl)
df3.insert(12,"log_std_v",fitstdv)
df3.insert(13,"log_rho",fitrho)
df3.insert(14,"fit_goodness",fitgoodness)
df3.insert(15,"std_gfp_correct",np.sqrt(1-np.power(np.array(df3.log_rho),2))*np.array(df3.log_std_gfp))
df3.to_csv("gated.csv")
#
#
# plt.plot([0,5],[0,5],[0,5],[0,10])
#
#
# for index,row in df.iterrows():
#     try:
#         meta,data = fcsparser.parse("../FCS/"+row.filename, meta_data_only=False, reformat_meta=True)
#         data.columns=[x.strip().replace('-', '_') for x in data.columns]
#         data=data[(data["SSC_H"]>np.exp(2.5)) & (data["SSC_A"]>0)&(data["FSC_H"]>np.exp(1.5))&(data["GFP_H"]>0)&(data["SSC_A"]>data["SSC_H"])&(data["SSC_H"]>data["FSC_H"])]
#         ax1=plt.subplot(131)
#         ax1.hist2d(np.log(data["FSC_H"]),np.log(data["SSC_H"]),bins=100,range=[[0,5],[0,5]]);
#         plt.plot([0,5],[0,5],[0,5],[0,10])
#         ax2=plt.subplot(132)
#         ax2.hist2d(np.log(data["SSC_A"]),np.log(data["SSC_H"]),bins=100,range=[[0,5],[0,5]]);
#         ax2.plot([0,5],[0,5])
#         ax3=plt.subplot(133)
#         ax3.hist2d(np.log(data["GFP_H"]),np.log(data["FSC_H"]),bins=100,range=[[0,5],[0,5]]);plt.show()
#         plt.suptitle(index)
#         print(row)
#         #data.plot(x=np.log("FSC_H"),y=np.log("SSC_H"),kind="hexbin",loglog=False)#,xlim=[1e-3,100],ylim=[1e-3,100]
#     except f:
#         print("fail at"+str(index))
