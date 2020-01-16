import pandas
%matplotlib inline
df=pandas.read_csv("volume_decomposed.csv")
df.iterrows

import fcsparser
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import warnings
warnings.simplefilter("default")
iters=df[df.fit_goodness>0.98].sort_values(by="log_rho",ascending=False).head(20).iterrows()
# for index, row in iters:
#     print("test")#shameless cherrypicking for example
#
#
#     data_dir='../FCS/'
#     example_file=row.filename#"Huseyin2019-09-26.0211.fcs"#"Huseyin2019-09-24.0072.fcs"#
#     meta, data = fcsparser.parse(data_dir + example_file, meta_data_only=False, reformat_meta=True)
#     data.columns=[x.strip().replace('-', '_') for x in data.columns]
#
#     df2=row#df[df.filename==example_file].iloc[0]
#     data["volume_decomposed"]=np.exp((np.log(data["GFP_H"])*df2["log_std_v"]-df2["log_std_gfp"]*np.log(data["FSC_H"])*df2["log_rho"]+df2["log_std_gfp"]*df2["log_mean_v_mean"]*df2["log_rho"])/df2["log_std_v"])
#     f, [ax,ax1] = plt.subplots(ncols=2,figsize=(10, 5),sharey=True)
#
#     f.suptitle(row.filename)
#     #ax.set_ylim([-1,5])
#     #ax1.set_ylim([-1,5])
#     #ax.set(xscale="log", yscale="log")
#     # Basic 2D density plot
#     sns.set_style("white")
#     sns.kdeplot(np.log(data.head(1000).FSC_H), np.log(data.head(1000).GFP_H),ax=ax)
#     sns.kdeplot(np.log(data.head(1000).FSC_H), np.log(data.head(1000).GFP_H), cmap="Blues", shade=True, shade_lowest=True,ax=ax )
#     sns.kdeplot(np.log(data.head(1000).FSC_H), np.log(data.head(1000).volume_decomposed),ax=ax1)
#     sns.kdeplot(np.log(data.head(1000).FSC_H), np.log(data.head(1000).volume_decomposed), cmap="Blues", shade=True, shade_lowest=True,ax=ax1 )
#     plt.show()
# "Huseyin2019-10-03.0125.fcs"
# "Huseyin2019-09-24.0234.fcs"
#df[df.filename=="Huseyin2019-09-24.0234.fcs"]
isinstance

iters2=df[df.filename.isin(["Huseyin2019-10-02.0024.fcs"])].sort_values(by="log_rho",ascending=False).head(10).iterrows()


def multivariate_gaussian(pos, mu, Sigma):
    """Return the multivariate Gaussian distribution on array pos.

    pos is an array constructed by packing the meshed arrays of variables
    x_1, x_2, x_3, ..., x_k into its _last_ dimension.

    """

    n = mu.shape[0]
    Sigma_det = np.linalg.det(Sigma)
    Sigma_inv = np.linalg.inv(Sigma)
    N = np.sqrt((2*np.pi)**n * Sigma_det)
    # This einsum call calculates (x-mu)T.Sigma-1.(x-mu) in a vectorized
    # way across all the input variables.
    fac = np.einsum('...k,kl,...l->...', pos-mu, Sigma_inv, pos-mu)

    return np.exp(-fac / 2) / N
for index, row in iters2:
    print("test")#shameless cherrypicking for example


    data_dir='../FCS/'
    example_file=row.filename#"Huseyin2019-09-26.0211.fcs"#"Huseyin2019-09-24.0072.fcs"#
    meta, data = fcsparser.parse(data_dir + example_file, meta_data_only=False, reformat_meta=True)
    data.columns=[x.strip().replace('-', '_') for x in data.columns]

    df2=row#df[df.filename==example_file].iloc[0]
    data["GFP_Decomposed"]=np.exp((np.log(data["GFP_H"])*df2["log_std_v"]-df2["log_std_gfp"]*np.log(data["FSC_H"])*df2["log_rho"]+df2["log_std_gfp"]*df2["log_mean_v_mean"]*df2["log_rho"])/df2["log_std_v"])
    f, [[ax,ax1],[ax2,ax3]] = plt.subplots(ncols=2,nrows=2,figsize=(10, 10),sharey=True,sharex=True,gridspec_kw = {'wspace':0, 'hspace':0})

    f.suptitle(row.filename)
    ax.set_ylim([2,4.5])
    ax.set_xlim([1,4])
    ax2.axvline(df2["log_mean_v_mean"])
    ax2.text(df2["log_mean_v_mean"]+0.02,1.5+1+0.15,'Context average',rotation=90)
    ax.text(df2["log_mean_v_mean"]+0.02,1.5+1+0.15,'Context average',rotation=90)
    ax.axvline(df2["log_mean_v_mean"])

    #ax1.set_ylim([-1,5])
    #ax.set(xscale="log", yscale="log")
    # Basic 2D density plot
    sns.set_style("white")
    data=data[(data["FSC_H"]>0)&(data["GFP_H"]>0)].head(100000)
    sns.kdeplot(np.log(data.FSC_H), np.log(data.GFP_H),ax=ax)
    sns.kdeplot(np.log(data.FSC_H), np.log(data.GFP_H), cmap="Blues", shade=True, shade_lowest=True,ax=ax )
    sns.kdeplot(np.log(data.FSC_H), np.log(data.GFP_Decomposed),ax=ax1)
    sns.kdeplot(np.log(data.FSC_H), np.log(data.GFP_Decomposed), cmap="Blues", shade=True, shade_lowest=True,ax=ax1 )
    N = 60
    X = np.linspace(1, 4, N)
    Y = np.linspace(2, 4.5, N)
    X, Y = np.meshgrid(X, Y)
    pos = np.empty(X.shape + (2,))
    pos[:, :, 0] = X
    pos[:, :, 1] = Y
    mu = np.array([row.log_mean_v, row.log_mean_gfp])
    Sigma = np.array([[ row.log_std_v**2 , row.log_rho*row.log_std_v*row.log_std_gfp], [row.log_rho*row.log_std_v*row.log_std_gfp,  row.log_std_gfp**2]])
    Z = multivariate_gaussian(pos, mu, Sigma)
    cset = ax2.contourf(X, Y, Z,cmap="Blues")
    test=np.sqrt(1-np.power(np.array(row.log_rho),2))*np.array(row.log_std_gfp)
    newmean=(df2["log_mean_gfp"]*df2["log_std_v"]-df2["log_std_gfp"]*df2["log_mean_v"]*df2["log_rho"]+df2["log_std_gfp"]*df2["log_mean_v_mean"]*df2["log_rho"])/df2["log_std_v"]
    mu = np.array([row.log_mean_v, newmean])
    Sigma = np.array([[ row.log_std_v**2 , 0], [0,  test**2]])
    Z = multivariate_gaussian(pos, mu, Sigma)
    ax3.contourf(X, Y, Z,cmap="Blues")
    ax3.axhline(newmean)
    ax3.text(1.2,newmean+.1,'new average')
    ax2.axhline(row.log_mean_gfp)
    ax2.text(1.2,row.log_mean_gfp+.1,'average')
    newmean
    row.log_mean_gfp
    ax2.set_ylabel("GFP_H")
    ax2.set_xlabel('FSC_H')
    ax3.set_xlabel('FSC_H')
    plt.savefig(row.filename+".svg")
    plt.show()
#
# fig = plt.figure()
# ax = fig.gca(projection='3d')       cmap=cm.viridis)
# Z = multivariate_gaussian(pos, mu, Sigma)
# cset = plt.contourf(X, Y, Z, zdir='z', offset=-0.15,)
#     # ax=data.plot.scatter(x='FSC_H',y='GFP_H')
#     # ax.set_yscale('log')
#     # ax.set_xscale('log')
#     # ax.set_ylim([0.1,1000])
#     # ax.set_title(row.filename)
#     #
#     # ax=data.plot.scatter(x='FSC_H',y='volume_decomposed')
#     # ax.set_yscale('log')
#     # ax.set_xscale('log')
#     # ax.set_ylim([0.1,1000])
#
#
#
#
#
#
# [data["GFP_H"],data["FSC_H"]]
