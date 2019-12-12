import pandas
import matplotlib.pyplot as plt
import os
import numpy as np
import re
%matplotlib inline

normal_df=pandas.read_csv('log_normal_fitted2.csv',index_col=0)
#normal_df=pandas.read_csv('normal_df.csv',index_col=0)
normal_df.strain=normal_df.strain.replace('CC118λPir','CC118Lpir')
normal_df.strain=normal_df.strain.replace('DH5α','DH5alpha')
normal_df.strain=normal_df.strain.replace('KT 2440','KT2440')
normal_df=normal_df.rename(columns={"iptg_concentration": "iptg","filename":"filename2"})
normal_df.plasmid=normal_df.plasmid.str.replace("-","_")
normal_df.plasmid=normal_df.plasmid.str.capitalize()
normal_df.plasmid=normal_df.plasmid.replace('Phif_pone','Phif_p1')
normal_df.iptg=normal_df.iptg.replace(0,1)

normal_df
normal_df2=normal_df.copy()
normal_df2["meanlogV"]=0
groups = normal_df2.groupby(['backbone','strain',"plasmid"])
for name, group in groups:
    group["meanlogV"]=group.log_mean_v.mean()
compare=pandas.merge(normal_df,normal_df.groupby(['backbone','strain']).log_mean_v.mean().to_frame(), on=["strain","backbone"], how='outer')
compare["rrpu"]=(compare["log_mean_gfp"]*compare["log_std_v"]-compare["log_std_gfp"]*compare["log_mean_v_x"]*compare["log_rho"]+compare["log_std_gfp"]*compare["log_mean_v_y"]*compare["log_rho"])/compare["log_std_v"]

huseyin_median=pandas.read_csv('huseyin_median.csv')
huseyin_median=huseyin_median.drop("Unnamed: 0",axis=1)
huseyin_median.plasmid=huseyin_median.plasmid.str.capitalize()
huseyin_median.plasmid=huseyin_median.plasmid.replace('Bmr3r1_b3','Bm3r1_b3')
huseyin_median.iptg=huseyin_median.iptg.replace(0,1)
huseyin_median

compare2=pandas.merge(compare,huseyin_median, on=["strain","plasmid","backbone","iptg"], how='outer',validate="1:1")

compare2.to_csv("median_and_mean.csv")
groups = compare2.query("backbone!=['Empty','EmptyRow']").groupby(['backbone','strain'])

fig, axs = plt.subplots(22,4,figsize=(30,8*22))
axs=axs.ravel()
for name, group in groups:
    #fig, [ax1,ax2,ax3,ax4] = plt.subplots(1,4,figsize=(30,8))
    background2=np.exp(group.log_mean_gfp.min())
    background= np.exp(group.query('plasmid=="1201"').rrpu.mean())#-np.exp(background.rrpu.mean()))
    constituant= group.query('plasmid=="1717"')
    standard=group.query('plasmid=="1818"')
    group2=group.query("plasmid==['1201','1717','Amer_f1', 'Beti_e1','Bm3r1_b1', 'Bm3r1_b2', 'Psra_r1']")
    plasms=group.groupby('plasmid')
    for id,[plas,plasm] in enumerate(plasms):
        print(id)
        plotter=(np.exp(plasm.rrpu)-background)/(np.exp(constituant.rrpu.mean())-background)
        plotter2=(np.exp(plasm.log_mean_gfp)-background)/(np.exp(constituant.log_mean_gfp.mean())-background)
        #plotter2=np.log((np.exp(plasm.log_mean_gfp)-np.exp(background.log_mean_gfp.mean()))/(np.exp(constituant.log_mean_gfp.mean())-np.exp(background.log_mean_gfp.mean())))
        axs[id+1].plot(plasm.iptg, np.exp(plasm.log_mean_v_x), marker='o', linestyle='-', ms=12, label=plas)
        axs[id+2].plot(plasm.iptg, plotter2, marker='o', linestyle='-', ms=12, label=plas)
        axs[id+3].plot(plasm.iptg, plotter, marker='o', linestyle='-', ms=12, label=plas)
        axs[id+4].plot(plasm.iptg, plasm.median_yfp, marker='o', linestyle='-', ms=12, label=name)
        axs[id+1].set_xscale('log')
        axs[id+2].set_xscale('log')
        axs[id+3].set_xscale('log')
        axs[id+4].set_xscale('log')
        #axs[id+1].set_yscale('log')
        # axs[id+2].set_yscale('log')
        # axs[id+3].set_yscale('log')
        # axs[id+4].set_yscale('log')
        axs[id+1].set_xlabel('iptg')
        axs[id+1].set_xlabel('iptg')
        axs[id+1].set_xlabel('iptg')
        axs[id+1].set_ylabel('FSC')
        axs[id+2].set_ylabel('not volume conditioned rpu')
        axs[id+3].set_ylabel('volume conditioned rpu')
        axs[id+1].legend()
        axs[id+2].legend()
        axs[id+3].legend()
        axs[id+4].legend()
        axs[id+1].set_ylim([0.,50])
        # axs[id+2].set_ylim([0.01,50])
        # axs[id+3].set_ylim([0.01,50])
        # axs[id+4].set_ylim([0.01,50])
    fig.suptitle(name)
    name2 = re.sub('[^0-9a-zA-Z]+', '_', str(name))
    #plt.savefig(name2+".png")
plt.show()


compare.plasmid.unique()
fig, [ax1,ax2,ax3] = plt.subplots(1,3,figsize=(30,8))
#ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
for name, group in groups:

    # background= group.query('plasmid=="1201"')
    # constituant= group.query('plasmid=="1717"')
    # standard=group.query('plasmid=="18181"')
    group2=group.query('plasmid=="1201"')
    plotter=(group2.rrpu-group3.rrpu.mean())
    ax1.plot(group2.iptg, group2.log_mean_v_x, marker='o', linestyle='-', ms=12, label=name)
    ax2.plot(group2.iptg, group2.log_mean_gfp, marker='o', linestyle='-', ms=12, label=name)
    ax3.plot(group2.iptg, group2.rrpu, marker='o', linestyle='-', ms=12, label=name)
ax1.legend()
ax1.set_xscale('log')
ax2.legend()
ax2.set_xscale('log')
ax3.legend()
ax3.set_xscale('log')
plt.show()

compare.corr()
