import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import re
import scipy.stats as st
import multiprocessing
#import pysftp
import os

#with pysftp.Connection('ruudstoof.com', username='uploader', password='myPass') as sftp:
#        print(sftp.listdir("upload/Huseyin_Tas_cello_in_putida_flow_cytometry_non_extended"))
#        with sftp.cd('/Huseyin_Tas_cello_in_putida_flow_cytometry_non_extended'):           # temporarily chdir to allcode
#            #sftp.put('/pycode/filename')  	# upload file to allcode/pycode on remote
#            sftp.get('remote_file')




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
        fluormin = np.percentile(fl,2)
        fluormax = np.percentile(fl,98)
    volmin = np.percentile(vl,2)
    volmax = np.percentile(vl,98)
    xx, yy = np.mgrid[fluormin:fluormax:160j, volmin:volmax:160j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([fl, vl])
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    axs2[id].imshow(np.rot90(f),aspect="equal", cmap=plt.cm.gist_earth_r,extent=[fluormin, fluormax, volmin, volmax])
    return [xx,yy,f]


def fit2binormal(kdedat):
    #fits a 2D -normal distribution to 2D histogram of log transformed data
    [xx,yy,f]=kdedat
    xy_mesh = [xx,yy]
    guess_vals = [np.max(f), 7, 7, 1, 1, 1]
    #guess_vals = [np.max(f), 2, 3, 1, 1, .6]
    # perform the fit, making sure to flatten the noisy data for the fit routine
    fit_params, cov_mat = curve_fit(gaussian_2d, xy_mesh, np.ravel(f), p0=guess_vals, maxfev=100000,bounds=([0,-3,-3,0,0,0],[10,10,10,5,5,1]))
    # calculate fit parameter errors from covariance matrix
    fit_errors = np.sqrt(np.diag(cov_mat))
    # manually calculate R-squared goodness of fit
    fit_residual = f - gaussian_2d(xy_mesh, *fit_params).reshape(np.outer(xx[:,0],yy[0]).shape)
    fit_Rsquared = 1 - np.var(fit_residual)/np.var(f)
    return [fit_params,fit_Rsquared]


def calc_stuff(it):
    fit=fit2binormal(kdedata( processfile(it//12,it%12)))
    return fit[0]


#with pysftp.Connection('ruudstoof.com', username='uploader', password='myPass') as sftp:
#        print(sftp.listdir("upload/Huseyin_Tas_cello_in_putida_flow_cytometry_non_extended"))
#        with sftp.cd('/Huseyin_Tas_cello_in_putida_flow_cytometry_non_extended'):           # temporarily chdir to allcode
#            #sftp.put('/pycode/filename')  	# upload file to allcode/pycode on remote
#            sftp.get('remote_file')

mydir='../FlowRepository_FR-FCM-ZZQM_files'
if os.path.isdir(mydir):
    usedir=mydir
else:
    from tkinter import filedialog
    from tkinter import *
    root = Tk()
    root.withdraw()
    yourdir = filedialog.askdirectory(title = "Select Measurement folder")
    usedir=yourdir

ls
# mydir='../FlowRepository_FR-FCM-ZZQM_files'
os.chdir(usedir)
test=os.getcwd()
test
sys("pwd")
import inspect
src_file_path = inspect.getfile(lambda: None)
src_file_path
ls
os.path.dirname(os.path.abspath("__file__"))
dir_path = os.path.dirname(os.path.realpath(__file__))
this_py_file = os.path.realpath(__file__)
# os.path
# dotdot=os.path.realpath(usedir+"/..")
names=np.genfromtxt("Michael_Jahn_2016.csv", delimiter=',',dtype="str")
names


import fcsparser
%matplotlib inline

names[0,1:]
import sys
file
reload(sys)
sys.setdefaultencoding('utf8')
path=usedir+"/"+files[0]
usedir
names[0,1]
#dotdot+"/FlowRepository_FR-FCM-ZZQM_files/"+
names[0,1:]
plt.figure(figsize=(30,15))
for id,name in enumerate(names[0,1:]):
    meta, data = fcsparser.parse("FlowRepository_FR-FCM-ZZQM_files/"+name, meta_data_only=False, reformat_meta=True)
    plt.scatter(data["FS Log"],data['FL 1 Log'],.01)


names[-1,1:-2]
plt.figure(figsize=(30,15))
for id,name in enumerate(names[-1,1:-2]):
    meta, data = fcsparser.parse("FlowRepository_FR-FCM-ZZQM_files/"+name, meta_data_only=False, reformat_meta=True)
    plt.scatter(data["FS Log"],data['FL 1 Log'],.01)
    #,['r.', 'g.', 'b.', 'y.','r.', 'g.', 'b.', 'y.','r.', 'g.', 'b.', 'y.'][id])


Chan1="SS Log"
Chan2="FL 1 Log"


meta, data = fcsparser.parse("FlowRepository_FR-FCM-ZZQM_files/"+names[-1,1], meta_data_only=False, reformat_meta=True)


plt.plot(data["FS Log"],data['FL 1 Log'],"b.")
#np.genfromtxt(file,encoding='utf8')






fig, axs= plt.subplots(nrows=21, ncols=3, figsize=(6,20), dpi=150)
fig2, axs2= plt.subplots(nrows=21, ncols=3, figsize=(6,20), dpi=150)
axs=axs.ravel()
axs2=axs2.ravel()
fits=[]
for id,path in enumerate(files):
    path=usedir+"/"+path
    meta, data = fcsparser.parse(path, meta_data_only=False, reformat_meta=True)
    print(id)
    try:
        axs[id].plot(data[Chan1],data[Chan2],"b.")
        print(id)
        axs[id].set_xscale('log')
        axs[id].set_yscale('log')
        data=np.array([data[Chan1],data[Chan2]]).transpose()
        data=data[data[:,0]>1600]#deletes all negative measurements
        data=data[data[:,1]>100]
        fluor=np.log(data[:,0])   #log transformes data
        vol=np.log(data[:,1])
        fits.append(fit2binormal(kdedata([fluor,vol])))
    except:
        print("fail at id: "+str(id))


fig.savefig('test_fsc_controls.png', bbox_inches='tight')#, dpi=1500)
fig2.savefig('test2_fsc_controls.png', bbox_inches='tight')#, dpi=1500)

fits[:1]



np.savetxt("fits.csv",np.array(np.array(fits)[:,0]),delimiter=",")
#pool = multiprocessing.Pool()
#out1 = pool.map(calc_stuff, range(0, 12*len(files_chopped)))
2+2
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
