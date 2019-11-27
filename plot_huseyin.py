import pandas
import matplotlib.pyplot as plt
%matplotlib inline

df3=pandas.read_csv('log_normal fitted.csv')

df3
df3.hist(bins=50,figsize=(10,10));
df3.corr()
examplemeasurement=df3.query("strain=='KT 2440'")
# plt.plot(examplemeasurement.iptg_concentration,examplemeasurement.log_mean_gfp,)
# plt.xscale("log")
# examplemeasurement.plot()

#test=examplemeasurement.plot(x="iptg_concentration",y="log_mean_gfp",kind="scatter",logx=True,xlim=[5/1.5,2000*1.5],yerr="log_std_gfp");
for id,group in examplemeasurement.groupby('plasmid'):
    group.plot(x="iptg_concentration",y="log_mean_v",kind="scatter",logx=True,xlim=[5/1.5,2000*1.5],yerr="log_std_v",title=id);
