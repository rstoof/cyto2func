import re
import matplotlib.pyplot as plt
import pandas


df=pandas.read_csv('standardised.csv',index_col=0)
groups = df.query("backbone!=['Empty','EmptyRow']").groupby(['backbone','strain'])

fig, axs = plt.subplots(23,4,figsize=(15,3*23))
axs=axs.ravel()
for name, group in groups:
    plasms=group.groupby('plasmid')
    for id,[plas,plasm] in enumerate(plasms):
        axs[3*id+1-1].plot(plasm.iptg, np.exp(plasm.log_mean_v_x), marker='o', linestyle='-', ms=6, label=name)
        axs[3*id+2-1].plot(plasm.iptg, plasm.volume_decomposed_log_mean_gfp, marker='o', linestyle='-', ms=6, label=name)
        axs[3*id+3-1].plot(plasm.iptg, plasm.rrpu, marker='o', linestyle='-', ms=6, label=name)

        axs[3*id+1-1].set_xscale('log')
        axs[3*id+2-1].set_xscale('log')
        axs[3*id+3-1].set_xscale('log')

        #axs[3*id+1-1].set_yscale('log')
        axs[3*id+2-1].set_yscale('log')
        axs[3*id+3-1].set_yscale('log')

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

        axs[3*id+1].set_ylim([0.001,50])
        axs[3*id+2].set_ylim([0.001,50])
        axs[3*id+3].set_ylim([0.001,50])
    fig.suptitle(name)
    name2 = re.sub('[^0-9a-zA-Z]+', '_', str(name))
    #plt.savefig(name2+".png")
plt.tight_layout()
plt.savefig("all_plasmids_png", dpi=300)
plt.show()
