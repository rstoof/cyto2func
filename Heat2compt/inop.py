x1 = np.log(iptgs[1:])
y1 = np.array(fittedgates[0]).transpose()[:,1:]
f = interp1d(x1, y1)


arr2=[]
y=np.linspace(-2,6,100);
for ipt in np.linspace(np.log(iptgs[1]),np.log(iptgs[-1]),100):
    pars=f(ipt);
    arr2.append(np.exp(-((y-pars[1])**2/(pars[3]**2)))/(2*np.pi*pars[3]));
fig = plt.figure(figsize=(8,8))
ax = fig.gca()
arr=np.array(arr2).transpose()

ax.imshow(arr, extent=[np.log(iptgs[1]), np.log(iptgs[-1]), -1, 5],origin="lower")
ax.clabel(cset, inline=0, fontsize=10)
ax.set_xlabel('Iptg')
ax.set_ylabel('Fluorescence')
ax.set_xlim( np.log(iptgs[1]), np.log(iptgs[-1]))
ax.set_ylim(-1, 5)

plt.title('2D Heatmap flour vs iptg autofluor')

plt.show()
x1 = np.log(iptgs[1:])
y1 = np.array(fittedgates[0+1]).transpose()[:,1:]
f = interp1d(x1, y1)


arr3=[]
y=np.linspace(-2,6,100);
for ipt in np.linspace(np.log(iptgs[1]),np.log(iptgs[-1]),100):
    pars=f(ipt);
    arr3.append(np.exp(-((y-pars[1])**2/(pars[3]**2)))/(2*np.pi*pars[3]));
fig = plt.figure(figsize=(8,8))
ax = fig.gca()
arr=np.array(arr3).transpose()

ax.imshow(arr, extent=[np.log(iptgs[1]), np.log(iptgs[-1]), -1, 5],origin="lower")
ax.clabel(cset, inline=0, fontsize=10)
ax.set_xlabel('Iptg')
ax.set_ylabel('Fluorescence')
ax.set_xlim( np.log(iptgs[1]), np.log(iptgs[-1]))
ax.set_ylim(-1, 5)

plt.title('2D Heatmap flour vs iptg standard')

plt.show()

x1 = np.log(iptgs[1:])
y1 = np.array(fittedgates[0+2]).transpose()[:,1:]
f = interp1d(x1, y1)


arr4=[]
y=np.linspace(-2,6,100);
for ipt in np.linspace(np.log(iptgs[1]),np.log(iptgs[-1]),100):
    pars=f(ipt);
    arr4.append(np.exp(-((y-pars[1])**2/(pars[3]**2)))/(2*np.pi*pars[3]));
fig = plt.figure(figsize=(8,8))
ax = fig.gca()
arr=np.array(arr4).transpose()

ax.imshow(arr, extent=[np.log(iptgs[1]), np.log(iptgs[-1]), -1, 5],origin="lower")
ax.clabel(cset, inline=0, fontsize=10)
ax.set_xlabel('Iptg')
ax.set_ylabel('Fluorescence')
ax.set_xlim( np.log(iptgs[1]), np.log(iptgs[-1]))
ax.set_ylim(-1, 5)

plt.title('2D Heatmap flour vs iptg Input promoter')

plt.show()

gate=7
x1 = np.log(iptgs[1:])
y1 = np.array(fittedgates[gate+2]).transpose()[:,1:]
f = interp1d(x1, y1)


arr5=[]
y=np.linspace(-2,6,100);
for ipt in np.linspace(np.log(iptgs[1]),np.log(iptgs[-1]),100):
    pars=f(ipt);
    arr5.append(np.exp(-((y-pars[1])**2/(pars[3]**2)))/(2*np.pi*pars[3]));
fig = plt.figure(figsize=(8,8))
ax = fig.gca()
arr=np.array(arr5).transpose()

ax.imshow(arr, extent=[np.log(iptgs[1]), np.log(iptgs[-1]), -1, 5],origin="lower")
ax.clabel(cset, inline=0, fontsize=10)
ax.set_xlabel('Iptg')
ax.set_ylabel('Fluorescence')
ax.set_xlim( np.log(iptgs[1]), np.log(iptgs[-1]))
ax.set_ylim(-1, 5)
gatename=(files_chopped[gate+3][0]).partition('_')[0]
plt.title(gatename)

plt.show()
