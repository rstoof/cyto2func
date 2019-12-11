import pandas
import matplotlib.pyplot as plt
import os
import numpy as np
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
normal_df2

groups = normal_df2.groupby(['backbone','strain',"plasmid"])
for name, group in groups:
    group["meanlogV"]=group.log_mean_v.mean()

groups.to_frame()
groups.log_mean_v.mean().to_frame()
compare=pandas.merge(normal_df,normal_df.groupby(['backbone','strain']).log_mean_v.mean().to_frame(), on=["strain","backbone"], how='outer')

compare
compare["rrpu"]=(compare["log_mean_gfp"]*compare["log_std_v"]-compare["log_std_gfp"]*compare["log_mean_v_x"]*compare["log_rho"]+compare["log_std_gfp"]*compare["log_mean_v_y"]*compare["log_rho"])/compare["log_std_v"]
compare

(loggfpmean*logvstd - loggfpstd*logvmean*rho + loggfpstd*meanlogvmean*rho)/logvstd

test=groups.log_mean_v.mean().unstack()
test
test["pSeva221"]

groups.unstack()
normal_df2
normal_df2




group.rrpu



compare.plasmid.unique()


plasms.head()

plas

groups = compare.query("backbone!=['Empty','EmptyRow']").groupby(['backbone','strain'])

#ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
import re
np.exp(plasm.rrpu)
(np.exp(plasm.rrpu)-np.exp(background.rrpu.mean()))/(np.exp(constituant.rrpu.mean())-np.exp(background.rrpu.mean()))
for name, group in groups:
    fig, [ax1,ax2,ax3] = plt.subplots(1,3,figsize=(30,8))
    background2=np.exp(group.log_mean_gfp.min())#-np.exp(background.rrpu.mean()))
    background= group.query('plasmid=="1201"')
    constituant= group.query('plasmid=="1717"')
    standard=group.query('plasmid=="1818"')
    group2=group.query("plasmid==['1201','1717','Amer_f1', 'Beti_e1','Bm3r1_b1', 'Bm3r1_b2', 'Psra_r1']")
    plasms=group2.groupby('plasmid')
    for plas,plasm in plasms:
        plotter=(np.exp(plasm.rrpu)-background2)/(np.exp(constituant.rrpu.mean())-background2)
        plotter2=(np.exp(plasm.log_mean_gfp)-background2)/(np.exp(constituant.log_mean_gfp.mean())-background2)
        #plotter2=np.log((np.exp(plasm.log_mean_gfp)-np.exp(background.log_mean_gfp.mean()))/(np.exp(constituant.log_mean_gfp.mean())-np.exp(background.log_mean_gfp.mean())))
        ax1.plot(plasm.iptg, np.exp(plasm.log_mean_v_x), marker='o', linestyle='-', ms=12, label=plas)
        ax2.plot(plasm.iptg, plotter2, marker='o', linestyle='-', ms=12, label=plas)
        ax3.plot(plasm.iptg, plotter, marker='o', linestyle='-', ms=12, label=plas+'new')
        # ax4.plot(plasm.iptg, plasm.rrpu, marker='o', linestyle='-', ms=12, label=plas+'new')
        ax1.set_xscale('log')
        ax2.set_xscale('log')
        ax3.set_xscale('log')
        ax1.set_yscale('log')
        ax2.set_yscale('log')
        ax3.set_yscale('log')
        ax1.set_xlabel('iptg')
        ax1.set_xlabel('iptg')
        ax1.set_xlabel('iptg')
        ax1.set_ylabel('FSC')
        ax2.set_ylabel('not volume conditioned rpu')
        ax3.set_ylabel('volume conditioned rpu')
        ax1.legend()
        ax2.legend()
        ax3.legend()
        # ax4.set_xscale('log')
    fig.suptitle(name)
    name2 = re.sub('[^0-9a-zA-Z]+', '_', str(name))
    plt.savefig(name2+".png")
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
