import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import re
import os
import scipy.stats as st

def gaussian_2d(xy_mesh, amp, xc, yc, sigma_x, sigma_y,rho):

    # unpack 1D list into 2D x and y coords
    (x, y) = xy_mesh

    # make the 2D Gaussian matrix
    gauss = amp*np.exp(-((x-xc)**2/(sigma_x**2)-(2*rho*(x-xc)*(y-yc))/(sigma_x*sigma_y)+(y-yc)**2/(sigma_y**2)))/(2*np.pi*sigma_x*sigma_y*np.sqrt(1-rho**2))

    # flatten the 2D Gaussian down to 1D
    return np.ravel(gauss)

def fit2binormal(kdedat):
    #fits a 2D -normal distribution to 2D histogram of log transformed data
    [xx,yy,f]=kdedat
    xy_mesh = [xx,yy]

    guess_vals = [np.max(f), 2, 3, 1, 1, .6]

    # perform the fit, making sure to flatten the noisy data for the fit routine
    fit_params, cov_mat = curve_fit(gaussian_2d, xy_mesh, np.ravel(f),
    p0=guess_vals,maxfev=100000,bounds=([0,-3,-3,0,0,0],[10,10,10,5,5,1]))

    # calculate fit parameter errors from covariance matrix
    fit_errors = np.sqrt(np.diag(cov_mat))

    # manually calculate R-squared goodness of fit
    fit_residual = f - gaussian_2d(xy_mesh, *fit_params).reshape(np.outer(xx[:,0],yy[0]).shape)
    fit_Rsquared = 1 - np.var(fit_residual)/np.var(f)
    return [fit_params,fit_Rsquared]
    #print('Fit R-squared: ', fit_Rsquared, '\n')
    #print('Fit Amplitude: ', fit_params[0], '\u00b1', fit_errors[0])
    #print('Fit X-Center:  ', fit_params[1], '\u00b1', fit_errors[1])
    #print('Fit Y-Center:  ', fit_params[2], '\u00b1', fit_errors[2])
    #print('Fit X-Sigma:   ', fit_params[3], '\u00b1', fit_errors[3])
    #print('Fit Y-Sigma:   ', fit_params[4], '\u00b1', fit_errors[4])
    #print('Fit Covariance:', fit_params[5], '\u00b1', fit_errors[5])
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

    xx, yy = np.mgrid[fluormin:fluormax:5j, volmin:volmax:5j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([fl, vl])
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    return [xx,yy,f]

def processfile( gate , iptg ):
    data=np.genfromtxt(direct+files_chopped[gate][iptg], delimiter='',skip_header=1)[:,[1,0]]
    data=data[data[:,0]>0]#deletes all negative measurements
    data=data[data[:,1]>0]
    fluor=np.log(data[:,0])   #log transformes data
    vol=np.log(data[:,1])
    #vol=vol[~np.isnan(fluor)]
    #fluor=fluor[~np.isnan(fluor)]
    #fluor=fluor[~np.isnan(vol)]
    #vol=vol[~np.isnan(vol)]

    return [fluor,vol]
def sort_nicely( l ):#""" Sort the given list in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    l.sort( key=alphanum_key)
    return l

import os
import re
direct='/home/ruud/Documents/upload/Distribution_Fits/renormed_data/'
files = os.listdir(direct)
files_alph=sort_nicely(files)
#print(files_alph)
n=12 #12 measurements each
files_chopped = [files_alph[i * n:(i + 1) * n] for i in range((len(files_alph) + n - 1) // n )]

def calc_stuff(it):
    fit=fit2binormal(kdedata( processfile(it//12,it%12)))
    return fit[0]
    #fitgate.append(fit[0])

import multiprocessing
#len(files_chopped)
pool = multiprocessing.Pool()
out1 = pool.map(calc_stuff, range(0, 12*len(files_chopped)))
np.array(out1).tofile("test.txt",",")
