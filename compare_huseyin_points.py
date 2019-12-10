import pandas
import matplotlib.pyplot as plt
import os
%matplotlib inline
df3=pandas.read_csv('log_normal_fitted.csv')
huseyin_median=pandas.read_csv('huseyin_median.csv')
huseyin_median=huseyin_median.drop("Unnamed: 0",axis=1)
huseyin_median.plasmid=huseyin_median.plasmid.str.capitalize()
huseyin_median.plasmid=huseyin_median.plasmid.replace('Bmr3r1_b3','Bm3r1_b3')
huseyin_median.iptg=huseyin_median.iptg.replace(0,1)
huseyin_median
just_medians=pandas.read_csv('just_medians.csv',index_col=0)
just_medians.strain=just_medians.strain.replace('CC118λPir','CC118Lpir')
just_medians.strain=just_medians.strain.replace('DH5α','DH5alpha')
just_medians.strain=just_medians.strain.replace('KT 2440','KT2440')
just_medians=just_medians.rename(columns={"iptg_concentration": "iptg","filename":"filename2"})
just_medians.plasmid=just_medians.plasmid.str.replace("-","_")
just_medians.plasmid=just_medians.plasmid.str.capitalize()
just_medians.plasmid=just_medians.plasmid.replace('Phif_pone','Phif_p1')
just_medians.iptg=just_medians.iptg.replace(0,1)


just_medians
testset1=huseyin_median.plasmid.unique()
testset2=just_medians.plasmid.unique()
set(testset1).intersection(set(testset2))
set(testset1).symmetric_difference(set(testset2))
testset1

testset2

huseyin_median
just_medians
len(huseyin_median.iptg.unique())
just_medians.iptg.unique()

compare=pandas.concat([huseyin_median,just_medians],ignore_index=True,axis=1,sort=False)
compare=pandas.concat([just_medians,huseyin_median],axis=1,join="inner")
compare=pandas.merge(just_medians,huseyin_median, on=["strain","plasmid","backbone","iptg"], how='outer')
compare[compare.plasmid=='1201'].strain.value_counts()

1717
1717
compare.corr()
