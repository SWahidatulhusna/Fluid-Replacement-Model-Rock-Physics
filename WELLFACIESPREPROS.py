import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

logs = pd.read_csv('Well_Tauk_AfterProcessing.csv')
print(logs.columns)
print("min", logs.Depth.min())
print("max", logs.Depth.max())
print("avg", logs.DT.mean())

logs = logs[(logs.Depth >= 6400) & (logs.Depth <= 6800)]

sand_cutoff = 0.50
brine_sand = ((logs.Vshale <= sand_cutoff) & (logs.SW >= 0.9))
oil_sand = ((logs.Vshale <= sand_cutoff) & (logs.SW < 0.9))
shale = (logs.Vshale > sand_cutoff)

temp_lfc = np.zeros(np.shape(logs.Vshale))
temp_lfc[brine_sand.values] = 1    # LFC will be 1 when ssb (brine sand flag) is True
temp_lfc[oil_sand.values] = 2      # LFC will be 2 when sso (oil sand flag) is True
temp_lfc[shale.values] = 4         # LFC will be 4 when sh (shale flag) is True
logs['LFC'] = temp_lfc             # Copy the temporary log temp_lfc into the DataFrame with name `LFC`

logs.to_csv('welltaukfacies_lfc.csv',index=False) # save the data for use in Part 2
print(np.shape(logs.Vshale))

string = "brine sst={0}, oil sst={1}, shale={2}"
data = (np.count_nonzero(brine_sand),
        np.count_nonzero(oil_sand),
        np.count_nonzero(shale))
print(string.format(*data))
print("LFC min: {0}, LFC max: {1}".format(logs.LFC.min(), logs.LFC.max()))

import matplotlib.colors as colors
#      0=undef   1=bri  2=oil   3=gas 4=shale
ccc = ['#B3B3B3','blue','green','red','#996633',]
cmap_facies = colors.ListedColormap(ccc[0:len(ccc)], 'indexed')

ztop = 6400; zbot=6800
ll=logs[(logs.Depth>=ztop) & (logs.Depth<=zbot)]

cluster=np.repeat(np.expand_dims(ll['LFC'].values,1), 100, 1)

f, ax = plt.subplots(nrows=1, ncols=4, figsize=(8, 12))
ax[0].plot(ll.Vshale, ll.Depth, '-g', label='Vsh')
ax[0].plot(ll.SW, ll.Depth, '-b', label='Sw')
ax[1].plot(ll.PhieTotal, ll.Depth, '-k', label='phi')
ax[1].plot(ll.NPHI, ll.Depth, '-', color='0.5')
ax[2].plot(ll.Vp0Vs, ll.Depth, '-', color='0.5')
im=ax[3].imshow(cluster, interpolation='none', aspect='auto',cmap=cmap_facies,vmin=0,vmax=4)
ax[1].fill_betweenx(ll.Depth, ll.PhieTotal,ll.NPHI , where=(ll.PhieTotal >= ll.NPHI), color='orange')
ax[1].fill_betweenx(ll.Depth, ll.PhieTotal, ll.NPHI, where=(ll.PhieTotal <= ll.NPHI), color='blue')
y2 = ll.Vshale
y1 = y2 * 0 + 0.5
ax[0].fill_betweenx(ll.Depth, y1, y2, where=(y1 >= y2), color='gold', linewidth=0)
ax[0].fill_betweenx(ll.Depth, y1, y2, where=(y1 <= y2), color='lime', linewidth=0)
cbar = plt.colorbar(im, ax=ax[3])
cbar.set_label('0=undef,1=brine,2=oil,3=gas,4=shale')
cbar.set_ticks(range(0,4+1)); cbar.set_ticklabels(range(0,4+1))
cbar.set_label((12*' ').join(['undef', 'ss brine', 'ss oil', 'ss gas', 'shale']))
cbar.set_ticks(range(0,1)); cbar.set_ticklabels('')

for i in range(len(ax)-1):
    ax[i].set_ylim(ztop,zbot)
    ax[i].invert_yaxis()
    ax[i].grid()
    ax[i].locator_params(axis='x', nbins=4)
ax[0].legend(fontsize='small', loc='lower right')
ax[0].set_xlabel("Vcl & Sw"),    ax[0].set_xlim(-.1,1.1)
ax[1].set_xlabel("Density & Nphi")
ax[2].set_xlabel("Vp/Vs"),
ax[3].set_xlabel('LFC')
ax[1].set_yticklabels([]); ax[2].set_yticklabels([]); ax[3].set_yticklabels([]); ax[3].set_xticklabels([])

plt.show()