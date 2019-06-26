import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import re
import os
import scipy.stats as st
import multiprocessing


mydir='/home/ruud/Documents/upload/Huseyin_Tas_cello_in_putida_flow_cytometry_non_extended/'
if os.path.isdir(mydir):
    usedir=mydir
else:
    from tkinter import filedialog
    from tkinter import *
    root = Tk()
    root.withdraw()
    yourdir = filedialog.askdirectory(title = "Select Measurement folder")
    usedir=yourdir


def sort_nicely( l ):#<3 wmil@stackexchange """ Sort the given list in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    l.sort( key=alphanum_key)
    return l

def gaussian_2d(xy_mesh, amp, xc, yc, sigma_x, sigma_y,rho):
    # unpack 1D list into 2D x and y coords
    (x, y) = xy_mesh
    # make the 2D Gaussian matrix
    gauss = amp*np.exp(-((x-xc)**2/(sigma_x**2)-(2*rho*(x-xc)*(y-yc))/(sigma_x*sigma_y)+(y-yc)**2/(sigma_y**2)))/(2*np.pi*sigma_x*sigma_y*np.sqrt(1-rho**2))
    # flatten the 2D Gaussian down to 1D
    return np.ravel(gauss)


files = os.listdir(usedir)
files_alph=sort_nicely(files)
n=12 #Assumes 12 measurements each
files_chopped = [files_alph[i * n:(i + 1) * n] for i in range((len(files_alph) + n - 1) // n )]


iptgs=[]
names=[]
for files_set in files_chopped:
    names.append((files_set[0]).partition('_')[0])
for conds in files_chopped[0]:
    iptgs.append((conds).partition('_')[2].partition('.')[0])


def processfile( gate , iptg ):
    data=np.genfromtxt(direct+files_chopped[gate][iptg], delimiter=',',skip_header=1)[:,[-5,5]]
    data=data[data[:,0]>0]#deletes all negative measurements
    data=data[data[:,1]>0]
    fluor=np.log(data[:,0])   #log transformes data
    vol=np.log(data[:,1])
    return [fluor,vol]


def kdedata(dat,usestandardrange=False): #makes a 2D histogram of log transformed data
    [fl,vl]=dat
    if usestandardrange:
        fluormin =np.log(.001)
        fluormax =np.log(100)
    else:
        fluormin = np.percentile(fl,2)
        fluormax = np.percentile(fl,98)
    volmin = np.percentile(vl,2)
    volmax = np.percentile(vl,98)
    xx, yy = np.mgrid[fluormin:fluormax:160j, volmin:volmax:160j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([fl, vl])
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    return [xx,yy,f]


def fit2binormal(kdedat):
    #fits a 2D -normal distribution to 2D histogram of log transformed data
    [xx,yy,f]=kdedat
    xy_mesh = [xx,yy]
    guess_vals = [np.max(f), 2, 3, 1, 1, .6]
    # perform the fit, making sure to flatten the noisy data for the fit routine
    fit_params, cov_mat = curve_fit(gaussian_2d, xy_mesh, np.ravel(f), p0=guess_vals, maxfev=100000,bounds=([0,-3,-3,0,0,0],[10,10,10,5,5,1]))
    # calculate fit parameter errors from covariance matrix
    fit_errors = np.sqrt(np.diag(cov_mat))
    # manually calculate R-squared goodness of fit
    fit_residual = f - gaussian_2d(xy_mesh, *fit_params).reshape(np.outer(xx[:,0],yy[0]).shape)
    fit_Rsquared = 1 - np.var(fit_residual)/np.var(f)
    return [fit_params,fit_Rsquared]


pool = multiprocessing.Pool()
out1 = pool.map(calc_stuff, range(0, 12*len(files_chopped)))
fittedgates = [out1[i * n:(i + 1) * n] for i in range((len(out1) + n - 1) // n )]





#now we assume that the standard is the second fittedgate

stdrdmeanflu=np.array(fittedgates)[1,:,1].mean()
stdrdmeanV=np.array(fittedgates)[1,:,2].mean()
stdrdstdflu=np.array(fittedgates)[1,:,3].mean()
stdrdstdV=np.array(fittedgates)[1,:,4].mean()
stdrdrho=np.array(fittedgates)[1,:,5].mean()
for gate in range(2,len(names)):
    print("start "+names[gate])
    #f, (axs) = plt.subplots(2, 6, sharex=True, sharey=True)
    #print(axs)
    [fluormin,fluormax]=[100,-100]
    for iptg in range(12):
        #print(iptgs[iptg])
        if iptg==0:
            ax=plt.subplot(12,1,1)
        else:
            ax=plt.subplot(12,1,iptg+1,sharex=ax,sharey=ax)
        [fluor,vol]=processfile( gate , iptg )
        dat=np.array([fluor,vol]).transpose()
        for subdat in dat:
            stdrdcondV=(subdat[1]*stdrdstdflu*stdrdrho-stdrdmeanV*stdrdstdflu*stdrdrho+stdrdmeanflu*stdrdstdV)/stdrdstdV
            test=(np.exp(subdat[0])- np.exp(np.array(fittedgates)[0,:,1].mean()))/( np.exp(stdrdcondV)- np.exp(np.array(fittedgates)[0,:,1].mean()))
            if test>0:
                #print(stdrdcondV)
                subdat[0]=np.log(test)
            else:
                subdat=[1,1]
        dat=dat[~np.isnan(dat[:,1])]
        dat=dat[~np.isnan(dat[:,0])]
        np.savetxt("FSC_normalise_"+names[gate]+"_at_"+str(iptgs[iptg])+"milmol.txt",dat)
        [fluor2,vol2]=dat.transpose()
        try:
            [xx,yy,f]=kdedata([fluor2,vol2],True)
        except:
            print("epic fail at"+str(iptg)+"_"+str(gate))
        fluormin =np.log(.001)#min(fluormin, np.percentile(fluor2,10))
        print(fluormin)
        fluormax =np.log(100)#min(max(fluormax, np.percentile(fluor2,90)),np.log(100))
        print(fluormax)
        #fluormax =np.percentile(fluor2,98)
        volmin = np.percentile(vol2,2)
        volmax = np.percentile(vol2,98)
        xx=np.exp(xx)
        yy=np.exp(yy)
        [fluormin, fluormax, volmin, volmax]=np.exp(np.array([fluormin, fluormax, volmin, volmax]))
        ax.scatter(np.exp(fluor2[::100]),np.exp(vol2[::100]),c="white",alpha=None,edgecolors='black')
        cfset = ax.contourf(xx,yy, f, cmap='coolwarm',alpha=.8)
        cset = ax.contour(xx,yy, f, colors='k')
        ax.clabel(cset, inline=0, fontsize=10)
        ax.set_xlabel('Standard Promoters equivalent\n at same FSC')
        ax.set_ylabel('FSC')
        ax.set_xlim(fluormin,fluormax)
        ax.set_ylim(volmin, volmax)
        ax.title.set_text('FSC condition-standardised\n '+names[gate]+' @ '+str(iptgs[iptg])+' muM IPTG')
        ax.set_xscale("log")
        ax.set_yscale("log")
    print("done "+names[gate])
    plt.tight_layout()
    from datetime import datetime
    now = datetime.now()
    date=now.strftime("%Y_%m_%d_%H%M%S")
    plt.savefig('Ruud_distribution_Gate_'+str(names[gate])+'_logsacle'+date+'.png',dpi=200)
    plt.show()
