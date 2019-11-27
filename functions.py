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


def processfile( gate , iptg , fname = False):
    if fname:
        data=np.genfromtxt(usedir+"/"+fname, delimiter=',',skip_header=1)[:,[-5,5]]
    else:
        data=np.genfromtxt(usedir+files_chopped[gate][iptg], delimiter=',',skip_header=1)[:,[-5,5]]
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
    axs2[id].imshow(np.rot90(f),aspect="auto", cmap=plt.cm.gist_earth_r,extent=[fluormin, fluormax, volmin, volmax])
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


def calc_stuff(it):
    fit=fit2binormal(kdedata( processfile(it//12,it%12)))
    return fit[0]
