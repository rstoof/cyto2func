import csv
import numpy as np
import matplotlib.pyplot as plt
import re
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

minvalfsc=[]
minvalgfp=[]
datearr=[]
droparr=[]
for index,row in df.iterrows():
    try:
        meta, data = fcsparser.parse("../"+row.filename, meta_data_only=False, reformat_meta=True)
        data.columns=[x.strip().replace('-', '_') for x in data.columns]

        datetime_object = datetime.strptime(meta['$DATE']+" "+meta['$BTIM'], '%Y-%b-%d %H:%M:%S')
        datearr.append(datetime_object)
        minvalfsc.append(data.SSC_H.quantile(.05))
        minvalgfp.append(data.GFP_H.quantile(.05))

    except:
        print("fail at"+str(index))
        droparr.append(index)
df=df.drop(droparr)

#plt.plot(datearr,minvalfsc)
#df2=df2.drop(["lowfsc","lowgfp","real_time","lowssc"],axis=1)
df2=df
#df2
#df
#len(minvalfsc)
df2.insert(5,"lowfsc",minvalfsc)
df2.insert(6,"lowgfp",minvalgfp)
df2.insert(7,"real_time",datearr)


df2.boxplot(column=["lowfsc","lowgfp"],by=['date',"plasmid"],figsize=(120,4),rot=-90)


plt.plot(minvalgfp)

fstfile=df.filename[0]


#WHY DO WE SEE negative VANLUES???
#HACK-FIX


    data



filterdata=data.query("GFP_H>0")
data.sort_values(by="GFP_H")
data
df
data["GFP_A"]

plt.hist(data["FSC_H"],100,range=[-5,50]);
plt.hist(data["GFP_A"],100,range=[-5,15]);
plt.hist(data["FSC_A"],100,range=[-50,150]);

print(df.plasmid.value_counts())
background=df.query("strain=='KT 2440' and plasmid=='pSeva221::1201' ")
standard=df.query("strain=='KT 2440' and plasmid=='pSeva221::1717' ")
measurement=df.query("strain=='KT 2440' and plasmid=='pSeva221::1818' ")

fig = plt.figure()
ax = fig.add_subplot(111)

df.plasmid.value_counts()

for id,fname in measurement.iterrows():
    print(fname.filename)

meta
seaborn.set()
data.columns
strains=df.strain.unique()

quantiles=[]
for strain in strains:
    tempquantile=[]
    all=df.query("strain=='"+strain+"' and plasmid=='pSeva221::1201' ")
    for index,row in all.iterrows():
        try:
            meta, data = fcsparser.parse("../"+row.filename, meta_data_only=False, reformat_meta=True)
            data.columns=[x.strip().replace('-', '_') for x in data.columns]
            tempquantile.append(data['SSC_A'].quantile(.05))
        except:
            print("fail at  "+row.filename)
    quantiles.append(tempquantile)
quantiles
quantiles2=[]
for strain in strains:
    tempquantile=[]
    all=df.query("strain=='"+strain+"' and plasmid=='pSeva221::1201' ")
    for index,row in all.iterrows():
        try:
            meta, data = fcsparser.parse("../"+row.filename, meta_data_only=False, reformat_meta=True)
            data.columns=[x.strip().replace('-', '_') for x in data.columns]
            tempquantile.append(data['GFP_A'].quantile(.05))
        except:
            print("fail at  "+row.filename)
    quantiles2.append(tempquantile)
quantiles2

for channel in ["GFP_H","FSC_H","SSC_H"]:
    if channel:
        for index,row in standard.iterrows():
            meta, data = fcsparser.parse("../"+row.filename, meta_data_only=False, reformat_meta=True)
            data.columns=[x.strip().replace('-', '_') for x in data.columns]
            filterdata=data.query(channel+"<200 and "+channel+">0")
            ax =seaborn.distplot(np.log(filterdata[channel]), hist=False, kde=True,label=str(row.iptg_concentration))
        ax.set_title(row.strain+"|"+row.plasmid)
        ax.legend(title="Iptg")
        plt.show()
for channel in ["GFP_A","FSC_A","SSC_A"]:
    if channel:
        for index,row in standard.iterrows():
            meta, data = fcsparser.parse("../"+row.filename, meta_data_only=False, reformat_meta=True)
            data.columns=[x.strip().replace('-', '_') for x in data.columns]
            filterdata=data.query(channel+"<200 and "+channel+">-100")
            ax =seaborn.distplot(filterdata[channel], hist=False, kde=True,label=str(row.iptg_concentration))
        ax.set_title(row.strain+"|"+row.plasmid)
        ax.legend(title="Iptg")
        plt.show()


df.plasmid.unique()

channels=["SSC_H","GFP_H"]
plasmid="pSeva221::AmeR_F1"#
for strain in strains:
    print("strain=='"+strain+"' and plasmid=='"+plasmid+"'")
    background=df.query("strain=='"+strain+"' and plasmid=='"+plasmid+"'")
    for index,row in background.iterrows():
        meta, data = fcsparser.parse("../"+row.filename, meta_data_only=False, reformat_meta=True)
        data.columns=[x.strip().replace('-', '_') for x in data.columns]
        #filterdata=data.query("SSC_A<100 and SSC_A>-50")
        ax = plt.plot(data[channels[0]],data[channels[1]],"+",alpha=0.01,label=str(row.iptg_concentration)+" corr:"+str(round(data.corr(method ='pearson')[channels[0]][channels[1]] ,2)))
        plt.yscale("log")
        plt.xscale("log")
        plt.xlabel(channels[0])
        plt.ylabel(channels[1])
    plt.title(row.strain+"|"+row.plasmid)
    plt.legend(title="Iptg")
    plt.show()


meta

corr
channels=["FSC_H","GFP_H"]
for strain in strains:
    measurement=df.query("strain== '"+strain+"' and plasmid=='pSeva221::1818' ")
    for index,row in measurement.iterrows():
        meta, data = fcsparser.parse("../"+row.filename, meta_data_only=False, reformat_meta=True)
        data.columns=[x.strip().replace('-', '_') for x in data.columns]
        data=data.query(channels[0]+">0 and "+channels[1]+">0")
        corr=np.corrcoef(np.log(data[channels[0]]),np.log(data[channels[1]]))[0,1]
        ax = plt.plot(np.log(data[channels[0]]),np.log(data[channels[1]]),"+",alpha=0.01,label=str(row.iptg_concentration)+" corr:"+str(corr))#+str(round(data.corr(method ='pearson')[channels[0]][channels[1]] ,2)))
        #plt.yscale("log")
        #plt.xscale("log")
        plt.xlabel(channels[0])
        plt.ylabel(channels[1])
    plt.title(row.strain+"|"+row.plasmid)
    plt.legend(title="Iptg")
    plt.show()



for strain in strains:
    standard=df.query("strain== '"+strain+"' and plasmid=='pSeva221::1717' ")
    for index,row in standard.iterrows():
        meta, data = fcsparser.parse("../"+row.filename, meta_data_only=False, reformat_meta=True)
        data.columns=[x.strip().replace('-', '_') for x in data.columns]
        #filterdata=data.query("SSC_A<100 and SSC_A>-50")
        tempquantileGFP=[]
        tempquantileSSC=[]
        all=df.query("strain=='"+strain+"' and plasmid=='pSeva221::1201' ")
        for index,row in all.iterrows():
            try:
                meta, data = fcsparser.parse("../"+row.filename, meta_data_only=False, reformat_meta=True)
                data.columns=[x.strip().replace('-', '_') for x in data.columns]
                tempquantileGFP.append(data['GFP_A'].quantile(.05))
                tempquantileSSC.append(data['SSC_A'].quantile(.05))
            except:
                print("fail at  "+row.filename)
        correctgfp=-1*np.mean(tempquantileGFP)
        correctssc=-1*np.mean(tempquantileSSC)
        print(correctgfp,correctssc)
        ax = plt.plot(data["GFP_A"]+correctgfp,data["SSC_A"]+correctssc,"+",alpha=0.01)
        plt.yscale("log")
        plt.xscale("log")
        plt.xlabel("GFP_A")
        plt.ylabel("SSC_A")
    plt.title(row.strain+"|"+row.plasmid)
    #plt.legend(title="Iptg")
    plt.show()



for index,row in measurement.iterrows():
    meta, data = fcsparser.parse("../"+row.filename, meta_data_only=False, reformat_meta=True)
    data.columns=[x.strip().replace('-', '_') for x in data.columns]
    #filterdata=data.query("SSC_A<100 and SSC_A>-50")
    ax = plt.plot(data["FSC_A"],data["SSC_A"],"+",alpha=0.01)
    plt.yscale("log")
    plt.xscale("log")
plt.show()

for index,row in measurement.iterrows():
    meta, data = fcsparser.parse("../"+row.filename, meta_data_only=False, reformat_meta=True)
    data.columns=[x.strip().replace('-', '_') for x in data.columns]
    #filterdata=data.query("SSC_A<100 and SSC_A>-50")
    ax = plt.plot(data["SSC_A"],data["SSC_H"],"+",alpha=0.01)
    plt.yscale("log")
    plt.xscale("log")
plt.show()
#ax.set_title(row.strain+"|"+row.plasmid)
#ax.legend(title="Iptg")
for index,row in measurement.iterrows():
    meta, data = fcsparser.parse("../"+row.filename, meta_data_only=False, reformat_meta=True)
    data.columns=[x.strip().replace('-', '_') for x in data.columns]
    #filterdata=data.query("SSC_A<100 and SSC_A>-50")
    ax = plt.plot(data["GFP_A"],data["GFP_H"],"+",alpha=0.01)
    plt.yscale("log")
    plt.xscale("log")
plt.show()

for index,row in measurement.iterrows():
    meta, data = fcsparser.parse("../"+row.filename, meta_data_only=False, reformat_meta=True)
    data.columns=[x.strip().replace('-', '_') for x in data.columns]
    filterdata=data.query("FSC_A<50 and FSC_A>-10")
    ax =seaborn.distplot(filterdata["FSC_A"], hist=False, kde=True,label=str(row.iptg_concentration))
ax.set_title(row.strain+"|"+row.plasmid)
ax.legend(title="Iptg")
plt.show()







#fitting to binormal######################################################################################################################################
from functions import gaussian_2d,kdedata,fit2binormal

import numpy as np

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
    plt.imshow(np.rot90(f),aspect="auto", cmap=plt.cm.gist_earth_r,extent=[fluormin, fluormax, volmin, volmax])
    plt.show()
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
channels=["FSC_H","GFP_H"]
testfile=df.head(1).filename[0]
testmeta,testdata=fcsparser.parse("../"+testfile, meta_data_only=False, reformat_meta=True)
testdata.columns=[x.strip().replace('-', '_') for x in testdata.columns]
testdata=testdata.query(channels[0]+">0 and "+channels[1]+">0")
testdata2=kdedata([np.log(testdata["GFP_H"]),np.log(testdata["FSC_H"])])
fit2binormal(testdata2)

minvalfsc=[]
minvalgfp=[]
datearr=[]
droparr=[]
fitarr=[]
for index,row in df.iterrows():
    try:
        meta, data = fcsparser.parse("../"+row.filename, meta_data_only=False, reformat_meta=True)
        data.columns=[x.strip().replace('-', '_') for x in data.columns]
        datetime_object = datetime.strptime(meta['$DATE']+" "+meta['$BTIM'], '%Y-%b-%d %H:%M:%S')
        datearr.append(datetime_object)
        minvalfsc.append(data.SSC_H.quantile(.05))
        minvalgfp.append(data.GFP_H.quantile(.05))
        data=data.query(channels[0]+">0 and "+channels[1]+">0")
        data2=kdedata([np.log(data["GFP_H"]),np.log(data["FSC_H"])])
        fits=fit2binormal(data2)
        fitarr.append(fits)
        print(fits)
    except:
        print("fail at"+str(index))
        droparr.append(index)
df

del df3
type(fitarr)
df3=df.copy()
df3=df3.drop(droparr)
df3
df3.insert(5,"lowfsc",minvalfsc)
df3.insert(6,"lowgfp",minvalgfp)
df3.insert(7,"real_time",datearr)

fitmufl=[fit[1] for fit in fitarr]
fitmuv=[fit[2] for fit in fitarr]
fitstdfl=[fit[3] for fit in fitarr]
fitstdv=[fit[4] for fit in fitarr]
fitrho=[fit[5] for fit in fitarr]
fitgoodness=[fit[6] for fit in fitarr]
fitmufl[1]


df3.insert(8,"log_mean_gfp",fitmufl)
df3.insert(9,"log_mean_v",fitmuv)
df3.insert(10,"log_std_gfp",fitstdfl)
df3.insert(11,"log_std_v",fitstdv)
df3.insert(12,"log_rho",fitrho)
df3.insert(13,"fit_goodness",fitgoodness)
df3.insert(14,"std_gfp_correct",np.sqrt(1-np.power(np.array(df3.log_rho),2))*np.array(df3.log_std_gfp))
df3.to_csv("log_normal_fitted.csv",index=False)

examplemeasurement=df3.query("strain=='KT 2440' and plasmid=='pSeva221::1818' ")

measurement


df3=pandas.read_csv('log_normal_fitted.csv')


df3
df3.hist(bins=50,figsize=(10,10));
df3.corr()
examplemeasurement=df3.query("strain=='KT 2440'")
plt.plot(examplemeasurement.iptg_concentration,examplemeasurement.log_mean_gfp,)
plt.xscale("log")
examplemeasurement.plot()

test=examplemeasurement.plot(x="iptg_concentration",y="log_mean_gfp",kind="scatter",logx=True,xlim=[5/1.5,2000*1.5],yerr="log_std_gfp");
grouped=examplemeasurement.groupby('plasmid')
plt.figure
len(grouped)
df3

df3.strain.value_counts()


fig, axs = plt.subplots(9,8, figsize=(30, 30), facecolor='w', edgecolor='k')
axs = axs.ravel()
for id,[name,group] in enumerate(examplemeasurement.groupby('plasmid')):
    print(id,name)
    group.plot(ax=axs[id],x="iptg_concentration",y="log_mean_gfp",kind="scatter",logx=True,xlim=[5/1.5,2000*1.5],yerr="log_std_gfp",title=name,color='r');
    group.plot(ax=axs[id],x="iptg_concentration",y="log_mean_gfp",kind="scatter",logx=True,xlim=[5/1.5,2000*1.5],yerr="std_gfp_correct",title=name);




for id,[name,group] in enumerate(examplemeasurement.groupby('plasmid')):
    print(id,name)
    group.plot(ax=axs[id],x="iptg_concentration",y="log_mean_v",kind="scatter",logx=True,xlim=[5/1.5,2000*1.5],yerr="log_std_v",title=name);



examplemeasurement.log_mean_v.mean()

from Equation import Expression


examplemeasurement

examplemeasurement.log_mean_v.mean()

2+2
test=examplemeasurement.plot(x="iptg_concentration",y="log_mean_v",kind="scatter",logx=True,xlim=[5/1.5,2000*1.5],yerr="log_std_v",by="plasmid");

test.set_xscale("log")
fitarr[:,1]

#droparr
fit2
# filterdata=data.query("GFP_H>0 and GFP_A>0")

ax.hist2d(data["FSC_H"],data["FSC_A"],100);

ax.set_xlim([-2,100])



filterdata=data.query("GFP_H>0 and GFP_A>0")
fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist2d(np.log(filterdata["GFP_H"]),np.log(filterdata["GFP_A"]),100);#,range=[[-1,4],[-2,2]]);









####################################################

for index,row in measurement.iterrows(): #plot gfp for "measurement"
    meta, data = fcsparser.parse("../"+row.filename, meta_data_only=False, reformat_meta=True)
    data.columns=[x.strip().replace('-', '_') for x in data.columns]
    filterdata=data.query("GFP_A<30 and GFP_A>-10")
    ax =seaborn.distplot(filterdata["GFP_A"], hist=False, kde=True,label=str(row.iptg_concentration))
ax.set_title(row.strain+"|"+row.plasmid)
ax.legend(title="Iptg")
plt.show()


for index,row in measurement.iterrows():#   plot ssc for "measurement"
    meta, data = fcsparser.parse("../"+row.filename, meta_data_only=False, reformat_meta=True)
    data.columns=[x.strip().replace('-', '_') for x in data.columns]
    filterdata=data.query("SSC_A<100 and SSC_A>-50")
    ax =seaborn.distplot(filterdata["SSC_A"], hist=False, kde=True,label=str(row.iptg_concentration))
ax.set_title(row.strain+"|"+row.plasmid)
ax.legend(title="Iptg")
plt.show()
