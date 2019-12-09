import pandas
import matplotlib.pyplot as plt
import os
%matplotlib inline

df3=pandas.read_csv('log_normal_fitted.csv')

df3
df3.hist(bins=50,figsize=(10,10));
df3.corr()

strains=df3.strain.unique()
for strain in strains:
    examplemeasurement=df3[(df3['strain']==strain)]
    # plt.plot(examplemeasurement.iptg_concentration,examplemeasurement.log_mean_gfp,)
    # plt.xscale("log")
    # examplemeasurement.plot()

    #test=examplemeasurement.plot(x="iptg_concentration",y="log_mean_gfp",kind="scatter",logx=True,xlim=[5/1.5,2000*1.5],yerr="log_std_gfp");
    fig, axs = plt.subplots(9,8, figsize=(50, 50), facecolor='w', edgecolor='k')
    fig.suptitle(strain)
    axs = axs.ravel()
    for id,[name,group] in enumerate(examplemeasurement.groupby('plasmid')):
        print(id,name)
        group.plot(ax=axs[id],x="iptg_concentration",y="log_mean_gfp",kind="scatter",logx=True,xlim=[5/1.5,2000*1.5],yerr="log_std_gfp",title=name,color='r');
        group.plot(ax=axs[id],x="iptg_concentration",y="log_mean_gfp",kind="scatter",logx=True,xlim=[5/1.5,2000*1.5],yerr="std_gfp_correct",title=name);
    fig.savefig(strain+".png")
