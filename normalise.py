min_to_back=False

import pandas
import numpy as np

#%matplotlib inline

df=pandas.read_csv('volume_decomposed.csv',index_col=0)



if min_to_back: #In our dataset the fluorescence sometimes goes below the autofluorescence, for comparison sometimeswe take the minimum value of fluorescence as the background
    df=pandas.merge(df,df.groupby(['backbone','strain']).volume_decomposed_log_mean_gfp.min().to_frame(), on=["strain","backbone"], how='outer',suffixes=("","_min"))
else:
    df=pandas.merge(df,df.query('plasmid=="1201"').groupby(['backbone','strain']).volume_decomposed_log_mean_gfp.mean().to_frame(), on=["strain","backbone"], how='outer',suffixes=("","_min"))
df=pandas.merge(df,df.query('plasmid=="1717"').groupby(['backbone','strain']).volume_decomposed_log_mean_gfp.mean().to_frame(), on=["strain","backbone"], how='outer',suffixes=("","_standard"))
df["rrpu"]=((np.exp(df.rrpu)-np.exp(df.volume_decomposed_log_mean_gfp_min))/(np.exp(df.volume_decomposed_log_mean_gfp_standard)-np.exp(df.volume_decomposed_log_mean_gfp_min)))
compactdf=df.copy()

compactdf=compactdf.drop(columns=['filename',"filename2","date","lowfsc","lowgfp","real_time","log_mean_v_y","volume_decomposed_log_mean_gfp_min","volume_decomposed_log_mean_gfp_standard"])

compactdf.to_csv("standardised.csv")

























































df["rrpu"]=(compare["log_mean_gfp"]*compare["log_std_v"]-compare["log_std_gfp"]*compare["log_mean_v_x"]*compare["log_rho"]+compare["log_std_gfp"]*compare["log_mean_v_y"]*compare["log_rho"])/compare["log_std_v"]

fig, axs = plt.subplots(23,4,figsize=(15,3*23))
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
        #print(id)
        plotter=(np.exp(plasm.rrpu)-background)/(np.exp(constituant.rrpu.mean())-background)
        plasm.newstandard=plotter
        plotter2=(np.exp(plasm.log_mean_gfp)-background)/(np.exp(constituant.log_mean_gfp.mean())-background)
        #plotter2=np.log((np.exp(plasm.log_mean_gfp)-np.exp(background.log_mean_gfp.mean()))/(np.exp(constituant.log_mean_gfp.mean())-np.exp(background.log_mean_gfp.mean())))
        axs[3*id+1-1].plot(plasm.iptg, np.exp(plasm.log_mean_v_x), marker='o', linestyle='-', ms=6, label=name)
        axs[3*id+2-1].plot(plasm.iptg, plotter2, marker='o', linestyle='-', ms=6, label=name)
        axs[3*id+3-1].plot(plasm.iptg, plotter, marker='o', linestyle='-', ms=6, label=name)

        axs[3*id+1-1].set_xscale('log')
        axs[3*id+2-1].set_xscale('log')
        axs[3*id+3-1].set_xscale('log')

        #axs[3*id+1-1].set_yscale('log')
        # axs[3*id+2-1].set_yscale('log')
        # axs[3*id+3-1].set_yscale('log')

        axs[3*id+1-1].set_xlabel('iptg')
        axs[3*id+2-1].set_xlabel('iptg')
        axs[3*id+3-1].set_xlabel('iptg')

        axs[3*id+2-1].set_ylabel('not volume conditioned rpu')
        axs[3*id+3-1].set_ylabel('volume conditioned rpu')
        axs[3*id+1-1].legend()
        axs[3*id+2-1].legend()
        axs[3*id+3-1].legend()

        axs[3*id+1-1].set_ylim([0.,50])
        axs[3*id+1-1].set_title("FSC_for:"+plas)
        axs[3*id+2-1].set_title("my_RPU_for:"+plas)
        axs[3*id+3-1].set_title("my_RRPU_for:"+plas)

        # axs[id+2].set_ylim([0.01,50])
        # axs[id+3].set_ylim([0.01,50])

    fig.suptitle(name)
    name2 = re.sub('[^0-9a-zA-Z]+', '_', str(name))
    #plt.savefig(name2+".png")
plt.tight_layout()
plt.savefig("all_plasmids"+".png", dpi=300)
plt.show()

df.newstandard
plasm.newstandard

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
