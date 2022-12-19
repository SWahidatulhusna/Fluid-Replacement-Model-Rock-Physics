import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def frm(vp1, vs1, rho1, rho_f1, k_f1, rho_f2, k_f2, k0, phi):
    vp1 = 1/(vp1*10**-3)*0.3
    vs1 = 1/(vs1*10**-3)*0.3
    mu1 = rho1 * vs1 ** 2.
    k_s1 = rho1 * vp1 ** 2 - (4. / 3.) * mu1

    # The dry rock bulk modulus
    kdry = (k_s1 * ((phi * k0) / k_f1 + 1 - phi) - k0) / ((phi * k0) / k_f1 + (k_s1 / k0) - 1 - phi)

    # Now we can apply Gassmann to get the new values
    k_s2 = kdry + (1 - (kdry / k0)) ** 2 / ((phi / k_f2) + ((1 - phi) / k0) - (kdry / k0 ** 2))
    rho2 = rho1 - phi * rho_f1 + phi * rho_f2
    mu2 = mu1
    vp2 = ((k_s2 + (4/3) * mu2) / rho2)**0.5
    vs2 = (mu2 / rho2)**0.5

    return vp2 * 1000, vs2 * 1000, rho2, k_s2

def vrh(volumes,k,mu):
    f = np.array(volumes).T
    k = np.resize(np.array(k),np.shape(f))
    mu = np.resize(np.array(mu),np.shape(f))

    k_u = np.sum(f*k, axis=1)
    k_l = 1. / np.sum(f/k, axis=1)
    mu_u = np.sum(f*mu, axis=1)
    mu_l = 1. / np.sum(f/mu, axis=1)
    k0 = (k_u+k_l) / 2.
    mu0 = (mu_u+mu_l) / 2.
    return k_u, k_l, mu_u, mu_l, k0, mu0

logs = pd.read_csv('welltaukfacies_lfc.csv')

rho_qz=2.65;  k_qz=37;  mu_qz=44    # mineral properties, quartz (i.e., sands)
rho_sh=2.81;  k_sh=15;  mu_sh=5     # mineral properties, clay (i.e., shales)
rho_b=1.09;   k_b=2.8               # fluid properties, brine
rho_o=0.78;   k_o=0.94              # fluid properties, oil
rho_g=0.25;   k_g=0.06              # fluid properties, gas



# mineral mixture bulk and shear moduli, k0 and mu0
shale = logs.Vshale.values
sand = 1 - shale - logs.PhieTotal.values
shaleN = shale / (shale+sand)  # normalized shale and sand volumes
sandN = sand / (shale+sand)
k_u, k_l, mu_u, mu_l, k0, mu0 = vrh([shaleN, sandN], [k_sh, k_qz], [mu_sh, mu_qz])

# fluid mixture bulk modulus, using the same vrh function but capturing the Reuss average (second output)
water = logs.SW.values
hc = 1 - water
tmp, k_fl, tmp, tmp, tmp, tmp = vrh([water, hc], [k_b, k_o], [0, 0])

# fluid mixture density
rho_fl = water*rho_b + hc*rho_o

Pr=(0.5-((logs.SDT/logs.DT)**2))/(1-((logs.SDT/logs.DT)**2))
logs['PR'] = Pr
vpb, vsb, rhobrine, kb = frm(logs.DT, logs.SDT, logs.RhoB, rho_fl, k_fl, rho_b, k_b, k0, logs.PhieTotal)
vpo, vso, rhooil, ko = frm(logs.DT, logs.SDT, logs.RhoB, rho_fl, k_fl, rho_o, k_o, k0, logs.PhieTotal)
vpg, vsg, rhogas, kg = frm(logs.DT, logs.SDT, logs.RhoB, rho_fl, k_fl, rho_g, k_g, k0, logs.PhieTotal)

sand_cutoff = 0.5
brine_sand = ((logs.Vshale <= sand_cutoff) & (logs.SW >= 1))
oil_sand = ((logs.Vshale <= sand_cutoff) & (logs.SW < 1))
shale = (logs.Vshale > sand_cutoff)

logs['VP'] = (1/(logs.DT*10**-3)*0.3)*10**3
logs['VS'] = (1/(logs.SDT*10**-3)*0.3)*10**3
logs['RHO'] = logs.RhoB



logs['VP_FRMB'] = (1/(logs.DT*10**-3)*0.3)*10**3
logs['VS_FRMB'] = (1/(logs.SDT*10**-3)*0.3)*10**3
logs['RHO_FRMB'] = logs.RhoB
logs['VP_FRMB'][brine_sand|oil_sand] = vpb[brine_sand|oil_sand]
logs['VS_FRMB'][brine_sand|oil_sand] = vsb[brine_sand|oil_sand]
logs['RHO_FRMB'][brine_sand|oil_sand] = rhobrine[brine_sand|oil_sand]
logs['IP_FRMB'] = logs.VP_FRMB*logs.RHO_FRMB
logs['IP_FRMB'][shale] = logs.IP
logs['VPVS_FRMB'] = logs.VP_FRMB/logs.VS_FRMB

logs['VP_FRMO'] = (1/(logs.DT*10**-3)*0.3)*10**3
logs['VS_FRMO'] = (1/(logs.SDT*10**-3)*0.3)*10**3
logs['RHO_FRMO'] = logs.RhoB
logs['VP_FRMO'][brine_sand|oil_sand] = vpo[brine_sand|oil_sand]
logs['VS_FRMO'][brine_sand|oil_sand] = vso[brine_sand|oil_sand]
logs['RHO_FRMO'][brine_sand|oil_sand] = rhooil[brine_sand|oil_sand]
logs['IP_FRMO'] = logs.VP_FRMO*logs.RHO_FRMO
logs['IP_FRMO'][shale] = logs.IP
logs['VPVS_FRMO'] = logs.VP_FRMO/logs.VS_FRMO

logs['VP_FRMG'] = (1/(logs.DT*10**-3)*0.3)*10**3
logs['VS_FRMG'] = (1/(logs.SDT*10**-3)*0.3)*10**3
logs['RHO_FRMG'] = logs.RhoB
logs['VP_FRMG'][brine_sand|oil_sand] = vpg[brine_sand|oil_sand]
logs['VS_FRMG'][brine_sand|oil_sand] = vsg[brine_sand|oil_sand]
logs['RHO_FRMG'][brine_sand|oil_sand] = rhogas[brine_sand|oil_sand]
logs['IP_FRMG'] = logs.VP_FRMG*logs.RHO_FRMG
logs['IP_FRMG'][shale]= logs.IP
logs['VPVS_FRMG'] = logs.VP_FRMG/logs.VS_FRMG

temp_lfc_b = np.zeros(np.shape(logs.Vshale))
temp_lfc_b[brine_sand.values | oil_sand.values] = 1  # LFC is 1 when either brine_sand (brine sand flag) or oil_sand (oil) is True
temp_lfc_b[shale.values] = 4                # LFC 4=shale
logs['LFC_B'] = temp_lfc_b

temp_lfc_o = np.zeros(np.shape(logs.Vshale))
temp_lfc_o[brine_sand.values | oil_sand.values] = 2  # LFC is now 2 when there's sand (brine_sand or oil_sand is True)
temp_lfc_o[shale.values] = 4                # LFC 4=shale
logs['LFC_O'] = temp_lfc_o

temp_lfc_g = np.zeros(np.shape(logs.Vshale))
temp_lfc_g[brine_sand.values | oil_sand.values] = 3  # LFC 3=gas sand
temp_lfc_g[shale.values] = 4                # LFC 4=shale
logs['LFC_G'] = temp_lfc_g

Imp1 = logs.IP.values
Imp2 = logs.IP_FRMB.values
Imp3 = logs.IP_FRMG.values
Imp4 = logs.IP_FRMO.values
t = logs.Checkshot.values
def Sintetic_Seismogram(Imp,t):
    dt = 0.001  # sampleing interval

    AI_tdom = np.interp(x=t, xp = logs.Checkshot, fp = Imp)    #resampling
    # again Rc calulation but in reampled time domain
    Rc_tdom = []
    for i in range(len(AI_tdom)-1):
        Rc_tdom.append((AI_tdom[i+1]-AI_tdom[i])/(AI_tdom[i]+AI_tdom[i+1]))
        # define function of ricker wavelet
    Rc_tdom.append(Rc_tdom[-1])
    def ricker(f, length, dt):
        t0 = np.arange(-length/2, (length-dt)/2, dt)
        Y = (1.0 - 2.0*(np.pi**2)*(f**2)*(t0**2)) * np.exp(-(np.pi**2)*(f**2)*(t0**2))
        return t0, Y
    f=20            #wavelet frequency
    length=0.512    #Wavelet vector length
    dt=dt           # Sampling prefer to use smiliar to resampled AI
    t0, w = ricker(f, length, dt) # ricker wavelet
    synthetic = np.convolve(w, Rc_tdom, mode='same')
    return synthetic

logs['SS_Original'] = Sintetic_Seismogram(Imp1,t)
logs['SS_Brine'] = Sintetic_Seismogram(Imp2,t)
logs['SS_Gas'] = Sintetic_Seismogram(Imp3,t)
logs['SS_Oil'] = Sintetic_Seismogram(Imp4,t)

import matplotlib.colors as colors
#      0=undef   1=bri  2=oil   3=gas 4=shale
ccc = ['#B3B3B3','blue','green','red','#996633',]
cmap_facies = colors.ListedColormap(ccc[0:len(ccc)], 'indexed')

ztop = 6750; zbot = 6800
ll = logs[(logs.Depth >= ztop) & (logs.Depth <= zbot)]

cluster=np.repeat(np.expand_dims(ll['LFC'].values,1),100,1)

f, ax = plt.subplots(nrows=1, ncols=7, figsize=(8, 25))

f.suptitle('Fluid Replacement Model : WELL-TAUK-1')

ax[0].plot(ll.Vshale, ll.Depth, '-g', label='Vsh')
ax[0].plot(ll.SW, ll.Depth, '-b', label='Sw')
ax[0].plot(ll.NPHI, ll.Depth, '-k', label='phi')

ax[1].plot(ll.VP_FRMG, ll.Depth, '-r',label='Gas')
ax[1].plot(ll.VP_FRMB, ll.Depth, '-b',label='brine')
ax[1].plot(ll.VP_FRMO, ll.Depth, '-k',label='Oil')
ax[1].plot(ll.VP, ll.Depth, '-', color='0.5',label='Original')

ax[2].plot(ll.VS_FRMG, ll.Depth, '-r',label='Gas')
ax[2].plot(ll.VS_FRMB, ll.Depth, '-b',label='brine')
ax[2].plot(ll.VS_FRMO, ll.Depth, '-k',label='Oil')
ax[2].plot(ll.VS, ll.Depth, '-', color='0.5',label='Original')

ax[3].plot(ll.RHO_FRMB, ll.Depth, '-b',label='brine')
ax[3].plot(ll.RHO_FRMG, ll.Depth, '-r',label='Gas')
ax[3].plot(ll.RHO_FRMO, ll.Depth, '-k',label='Oil')
ax[3].plot(ll.RHO,ll.Depth, '-',color='0.5',label='Original')

ax[4].plot(ll.IP_FRMG, ll.Depth, '-r',label='Gas')
ax[4].plot(ll.IP_FRMB, ll.Depth, '-b',label='brine')
ax[4].plot(ll.IP_FRMO, ll.Depth, '-k',label='Oil')
ax[4].plot(ll.IP, ll.Depth, '-', color='0.5',label='Original')

ax[5].plot(ll.VPVS_FRMG, ll.Depth, '-r',label='Gas')
ax[5].plot(ll.VPVS_FRMB, ll.Depth, '-b',label='brine')
ax[5].plot(ll.VPVS_FRMO, ll.Depth, '-k',label='Oil')
ax[5].plot(ll.Vp0Vs, ll.Depth, '-', color='0.5',label='Original')

ax[6].plot(ll.SS_Original,ll.Depth, '-',label='Original')
ax[6].fill_betweenx(ll.Depth,ll.SS_Original, 0,  ll.SS_Original > 0,  color='0.5')
ax[6].plot(ll.SS_Brine, ll.Depth, '-b',label='brine')
ax[6].fill_betweenx(ll.Depth,ll.SS_Brine, 0,  ll.SS_Brine > 0,  color='b')
ax[6].plot(ll.SS_Gas, ll.Depth, '-r',label='Gas')
ax[6].fill_betweenx(ll.Depth,ll.SS_Gas, 0,  ll.SS_Gas > 0,  color='r',)
ax[6].plot(ll.SS_Oil, ll.Depth, '-k',label='Oil')
ax[6].fill_betweenx(ll.Depth,ll.SS_Oil, 0,  ll.SS_Oil > 0,  color='k')

for i in range(7):
    ax[i].set_ylim(ztop,zbot)
    ax[i].invert_yaxis()
    ax[i].grid()
    ax[i].locator_params(axis='x', nbins=4)
ax[0].legend(fontsize='small', loc='lower right')
ax[0].set_xlabel("Vcl/phi/Sw")
ax[1].legend(fontsize='small', loc='lower right')
ax[1].set_xlabel("Vp [m/s]")
ax[2].legend(fontsize='small', loc='lower right')
ax[2].set_xlabel("Vs[m/s]")
ax[3].legend(fontsize='small', loc='lower right')
ax[3].set_xlabel('Rho [gr/cc]')
ax[4].legend(fontsize='small', loc='lower right')
ax[4].set_xlabel("Ip [m/s*g/cc]")
ax[5].legend(fontsize='small', loc='lower right')
ax[5].set_xlabel("Vp/Vs")
ax[6].legend(fontsize='small', loc='lower right')
ax[6].set_xlabel('Sintetik Seismogram')
ax[6].legend(fontsize='small', loc='lower right')
ax[6].set_xlabel('Sintetik Seismogram FRM')
ax[6].set_xlabel('Sintetik Seismogram Brine')
ax[6].set_xlabel('Sintetik Seismogram Gas')
ax[6].set_xlabel('Sintetik Seismogram Oil')
ax[1].set_yticklabels([])
ax[2].set_yticklabels([])
ax[3].set_yticklabels([])
ax[4].set_yticklabels([])
ax[5].set_yticklabels([])
ax[6].set_yticklabels([])



f, ax = plt.subplots(nrows=1, ncols=4, sharey=True, sharex=True, figsize=(16, 4))

ax[0].scatter(ll.IP,ll.Vp0Vs,20,ll.LFC,marker='o',edgecolors='none',alpha=0.5,cmap=cmap_facies,vmin=0,vmax=4)
ax[1].scatter(ll.IP_FRMB,ll.VPVS_FRMB,20,ll.LFC_B,marker='o',edgecolors='none',alpha=0.5,cmap=cmap_facies,vmin=0,vmax=4)
ax[2].scatter(ll.IP_FRMO,ll.VPVS_FRMO,20,ll.LFC_O,marker='o',edgecolors='none',alpha=0.5,cmap=cmap_facies,vmin=0,vmax=4)
ax[3].scatter(ll.IP_FRMG,ll.VPVS_FRMG,20,ll.LFC_G,marker='o',edgecolors='none',alpha=0.5,cmap=cmap_facies,vmin=0,vmax=4)
ax[0].set_xlabel("Ip [m/s*g/cc]")
ax[0].set_ylabel("Vp/Vs")
ax[0].set_title('original data')

ax[1].set_xlabel("Ip [m/s*g/cc]")
ax[1].set_ylabel("Vp/Vs")
ax[1].set_title('FRM to brine')

ax[2].set_xlabel("Ip [m/s*g/cc]")
ax[2].set_ylabel("Vp/Vs")
ax[2].set_title('FRM to oil')

ax[3].set_xlabel("Ip [m/s*g/cc]")
ax[3].set_ylabel("Vp/Vs")
ax[3].set_title('FRM to gas')

for i in ax: i.grid()
plt.show()