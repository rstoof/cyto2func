import pandas

normal_df=pandas.read_csv('log_normal_fitted2.csv',index_col=0)
gated=pandas.read_csv('gated.csv',index_col=0)
gated.columns
merged=pandas.merge(normal_df,gated,on=['filename', 'date', 'strain', 'plasmid', 'backbone','iptg_concentration','real_time'])
merged.corr()


import pandas
import matplotlib.pyplot as plt
import os
import numpy as np
import re
%matplotlib inline

normal_df=pandas.read_csv('gated.csv',index_col=0)
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

#compare3.query('plasmid=="1717"').groupby(['backbone','strain']).log_mean_gfp
#compare2["newstandard"]=0


cheeky=True
if cheeky:
    compare3=pandas.merge(compare2,compare2.groupby(['backbone','strain']).rrpu.min().to_frame(), on=["strain","backbone"], how='outer',suffixes=("","_min"))
else:
    compare3=pandas.merge(compare2,compare2.query('plasmid=="1201"').groupby(['backbone','strain']).rrpu.mean().to_frame(), on=["strain","backbone"], how='outer',suffixes=("","_min"))
compare4=pandas.merge(compare3,compare2.query('plasmid=="1717"').groupby(['backbone','strain']).rrpu.mean().to_frame(), on=["strain","backbone"], how='outer',suffixes=("","_standard"))
compare4["newstandard"]=((np.exp(compare4.rrpu)-np.exp(compare4.rrpu_min))/(np.exp(compare4.rrpu_standard)-np.exp(compare4.rrpu_min)))
compare5=compare4.copy()

compare5.columns
compare5=compare5.drop(columns=['filename',"filename2","date","lowfsc","lowgfp","real_time","log_mean_v_y","rrpu_min","rrpu_standard"])

compare5.to_csv("standardised_rrpu2"+("_cheeky" if cheeky else "")+"_gated"+".csv")


groups = compare5.query("backbone!=['Empty','EmptyRow']").groupby(['backbone','strain'])
#compare5[['rrpu','rrpu_min','rrpu_standard']]

# compare5.newstandard
# compare5.newstandard.hist(bins=30)
fig, axs = plt.subplots(23,4,figsize=(15,4*23))
axs=axs.ravel()
for name, group in groups:
    #fig, [ax1,ax2,ax3,ax4] = plt.subplots(1,4,figsize=(30,8))

    plasms=group.groupby('plasmid')
    for id,[plas,plasm] in enumerate(plasms):
        axs[4*id+1-1].plot(plasm.iptg, np.exp(plasm.log_mean_v_x), marker='o', linestyle='-', ms=6, label=name)
        axs[4*id+2-1].plot(plasm.iptg, plasm.rrpu, marker='o', linestyle='-', ms=6, label=name)
        axs[4*id+3-1].plot(plasm.iptg, plasm.newstandard, marker='o', linestyle='-', ms=6, label=name)
        axs[4*id+4-1].plot(plasm.iptg, plasm.median_yfp, marker='o', linestyle='-', ms=6, label=name)
        axs[4*id+1-1].set_xscale('log')
        axs[4*id+2-1].set_xscale('log')
        axs[4*id+3-1].set_xscale('log')
        axs[4*id+4-1].set_xscale('log')
        #axs[4*id+1-1].set_yscale('log')
        axs[4*id+2-1].set_yscale('log')
        axs[4*id+3-1].set_yscale('log')
        axs[4*id+4-1].set_yscale('log')
        axs[4*id+1-1].set_xlabel('iptg')
        axs[4*id+2-1].set_xlabel('iptg')
        axs[4*id+3-1].set_xlabel('iptg')
        axs[4*id+4-1].set_ylabel('FSC')
        axs[4*id+2-1].set_ylabel('not volume conditioned rpu')
        axs[4*id+3-1].set_ylabel('volume conditioned rpu')
        axs[4*id+1-1].legend()
        axs[4*id+2-1].legend()
        axs[4*id+3-1].legend()
        axs[4*id+4-1].legend()
        axs[4*id+1-1].set_ylim([0.,50])
        axs[4*id+1-1].set_title("FSC_for:"+plas)
        axs[4*id+2-1].set_title("my_RPU_for:"+plas)
        axs[4*id+3-1].set_title("my_RRPU_for:"+plas)
        axs[4*id+4-1].set_title("huseyin_RPU_for:"+plas)
        axs[4*id+1].set_ylim([0.001,50])
        axs[4*id+2].set_ylim([0.001,50])
        axs[4*id+3].set_ylim([0.001,50])
    fig.suptitle(name)
    name2 = re.sub('[^0-9a-zA-Z]+', '_', str(name))
    #plt.savefig(name2+".png")
plt.tight_layout()
plt.savefig("all_plasmids_"+str(cheeky)+"gated"+".png", dpi=300)
plt.show()
