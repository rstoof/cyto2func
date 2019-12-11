import pandas
import matplotlib.pyplot as plt
import os
import numpy as np
%matplotlib inline

huseyin_median=pandas.read_csv('huseyin_median.csv')
huseyin_median=huseyin_median.drop("Unnamed: 0",axis=1)
huseyin_median.plasmid=huseyin_median.plasmid.str.capitalize()
huseyin_median.plasmid=huseyin_median.plasmid.replace('Bmr3r1_b3','Bm3r1_b3')
huseyin_median.iptg=huseyin_median.iptg.replace(0,1)
huseyin_median
# just_medians=pandas.read_csv('just_medians.csv',index_col=0)
# just_medians.strain=just_medians.strain.replace('CC118λPir','CC118Lpir')
# just_medians.strain=just_medians.strain.replace('DH5α','DH5alpha')
# just_medians.strain=just_medians.strain.replace('KT 2440','KT2440')
# just_medians=just_medians.rename(columns={"iptg_concentration": "iptg","filename":"filename2"})
# just_medians.plasmid=just_medians.plasmid.str.replace("-","_")
# just_medians.plasmid=just_medians.plasmid.str.capitalize()
# just_medians.plasmid=just_medians.plasmid.replace('Phif_pone','Phif_p1')
# just_medians.iptg=just_medians.iptg.replace(0,1)
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
testset1=huseyin_median.plasmid.unique()
testset2=normal_df.plasmid.unique()
set(testset1).intersection(set(testset2))
set(testset1).symmetric_difference(set(testset2))
testset1

testset2

huseyin_median
normal_df
len(huseyin_median.iptg.unique())
normal_df.iptg.unique()

compare=pandas.concat([huseyin_median,normal_df],ignore_index=True,axis=1,sort=False)
compare=pandas.concat([normal_df,huseyin_median],axis=1,join="inner")
compare=pandas.merge(normal_df,huseyin_median, on=["strain","plasmid","backbone","iptg"], how='outer',validate="1:1")
compare[compare.plasmid=='1201'].strain.value_counts()
compare

1717
1717
compare.corr()

compare.plot(x="median_yfp",y=np.exp("log_mean_gfp"), kind = 'scatter')

compare.plot(x="median_yfp",y="log_mean_gfp", kind = 'scatter',logx=True,xlim=[0.001,100])

groups = compare.groupby(['backbone','strain'])

# Plot
fig, ax = plt.subplots(figsize=(20,10))
ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
for name, group in groups:
    ax.plot(group.median_yfp, group.log_mean_gfp, marker='o', linestyle='', ms=12, label=name)
ax.legend()
ax.set_xscale('log')
plt.show()


compare
groups = compare.groupby(['backbone','strain',"iptg"])
for name, group in groups:



group.log_mean_gfp
groups = compare.query('plasmid=="1201"')

fig, [ax1,ax2,ax3] = plt.subplots(1,3,figsize=(30,8))
#ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
for name, group in groups:
    ax1.plot(group.iptg, group.log_mean_v, marker='o', linestyle='-', ms=12, label=name)
    ax2.plot(group.iptg, group.log_mean_gfp, marker='o', linestyle='-', ms=12, label=name)
ax1.legend()
ax1.set_xscale('log')
ax2.legend()
ax2.set_xscale('log')
plt.show()

group.log_mean_gfp
groups = compare.query('plasmid=="1717"')
groups = groups.groupby(['backbone','strain'])
fig, ax = plt.subplots(figsize=(20,10))
ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
for name, group in groups:
    ax.plot(group.iptg, group.log_mean_gfp, marker='o', linestyle='-', ms=12, label=name)
ax.legend()
ax.set_xscale('log')
plt.show()



group.log_mean_gfp
groups = compare.query('plasmid=="1818"')
groups = groups.groupby(['backbone','strain'])
fig, ax = plt.subplots(figsize=(20,10))
ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
for name, group in groups:
    ax.plot(group.iptg, group.log_mean_gfp, marker='o', linestyle='-', ms=12, label=name)
ax.legend()
ax.set_xscale('log')
plt.show()
groups = compare.groupby(['plasmid'])

# Plot
fig, ax = plt.subplots(figsize=(20,10))
ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
for name, group in groups:
    ax.plot(group.iptg, group.log_mean_gfp, marker='o', linestyle='', ms=12, label=name)
ax.legend()
ax.set_xscale('log')
plt.show()

compare.plot(x="median_yfp",y="log_mean_gfp", kind = 'scatter',logx=True,xlim=[0.001,100])

fig, ax = plt.subplots(figsize=(20,10))
ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
for name, group in groups:
    ax.hist(group.log_mean_gfp,label=name,bins=100)
ax.legend()

plt.show()


compare.fit_goodness.hist(bins=200)
plt.savefig('output.png')
plt.plot(compare.median_yfp,np.exp(compare.log_mean_gfp),".")
plt.xscale("log")
plt.yscale("log")
np.exp(10)
