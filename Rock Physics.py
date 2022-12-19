import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

file = './WELL_TAUK-1.las'
data = np.loadtxt(file,skiprows=38)


memories= ['Depth','Checkshot','RhoB','GR','NPHI','DT','PHIE_1','PHIT_1','SDT','SW']
data=pd.DataFrame(data,columns=memories)
data= data[['Depth','Checkshot','RhoB','GR','NPHI','DT','PHIE_1','PHIT_1','SDT','SW']]
data=data.values

rows,cols= 1,9
fig,ax= plt.subplots(nrows=rows,ncols=cols,figsize=(12,10),sharey=True)
nmemories = ['Checkshot','RhoB','GR','NPHI','DT','PHIE_1','PHIT_1','SDT','SW']
for i in range(cols):
    ax[i].plot(data[:,i+1],data[:,0])
    ax[i].set_ylim(max(data[:,0]),min(data[:,0]))
    ax[i].minorticks_on()
    ax[i].grid(which='major',linestyle='-',linewidth='0.5',color='green')
    ax[i].grid(which='minor',linestyle= ':',linewidth= '0.5',color='black' )
    ax[i].set_title('%s'%nmemories[i])
# y2=data[:,2]
# y1=y2*0+100
# ax[1].set_title('%s'%nmemories[1])
# ax[1].fill_betweenx(data[:,0],y1,y2,where=(y1>=y2),color='gold',linewidth=0)
# ax[1].fill_betweenx(data[:,0],y1,y2,where=(y1<=y2),color='lime',linewidth=0)
# plt.subplots_adjust(wspace=0)
# ax[2].plot(data[:,8], data[:,0])
# ax[2].plot(data[:,3], data[:,0])
# ax[2].set_ylim(max(data[:,0]), min(data[:,0]))
# ax[2].minorticks_on()
# ax[2].grid(which='major', linestyle='-', linewidth='0.5', color='green')
# ax[2].grid(which='minor', linestyle=':', linewidth='0.5', color='black')
# ax[2].set_title('Porositas Total & NPhi')
# ax[2].fill_betweenx(data[:,0], data[:,8], data[:,3], where=(data[:,8] >= data[:,3]), color='red')
# ax[2].fill_betweenx(data[:,0], data[:,8], data[:,3], where=(data[:,8] <= data[:,3]), color='blue')

vshale= (data[:,3]-min(data[:,3]))/(max(data[:,3]-min(data[:,3])))
print(vshale)

plt.show()
