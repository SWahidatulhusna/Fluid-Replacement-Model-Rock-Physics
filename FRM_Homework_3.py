# Spesifikasi: Fluid Replacement Modelling
# Versi      : HomeWork 3 (Gassman & Greenberg Castagna)
# Last Edited: 2020-10-19
# Programmer : Sabda Wahidatulhusna

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

def Vs_GreenbergCastagna(vp,rho,Km,mum,kfl,kbr,phi,rhowet):
    vp= 1/(vp*10**-3)*0.3
    M= (vp**2)*rho
    Mm= Km+(4/3)*mum
    d=M/(Mm-M)-kfl/(phi*(Mm-kfl))+kbr/(phi*(Mm-kbr))
    Mwet=d*(Mm/(1+d))
    Vpwet=(Mwet/rhowet)**0.5
    Ac = 0.80416
    Bc = -0.85588
    Vswet = Ac*Vpwet+Bc
    Vs = Vswet*(rhowet/rho)**0.5
    return Vs

logs = pd.read_csv('welltaukfacies_lfc.csv')

rho_qz=2.65;  k_qz=37;  mu_qz=44    # mineral properties, quartz (i.e., sands)
rho_b=1.1;   k_b=2.6               # fluid properties, brine
rho_g=0.1;   k_g=0.02              # fluid properties, gas


# mineral mixture bulk and shear moduli, k0 and mu0
shale = logs.Vshale.values
sand = 1 - shale - logs.NPHI.values
shaleN = shale / (shale+sand)  # normalized shale and sand volumes
sandN = sand / (shale+sand)
k_u, k_l, mu_u, mu_l, k0, mu0 = vrh([sandN], [k_qz], [mu_qz])

# fluid mixture bulk modulus, using the same vrh function but capturing the Reuss average (second output)
water = logs.SW.values
hc = 1 - logs.SW0.values
tmp, k_fl, tmp, tmp, tmp, tmp = vrh([water, hc], [k_b, k_g], [0, 0])

# fluid mixture density
rho_fl = water*rho_b + hc*rho_g
vpg, vsg, rhogas, kg = frm(logs.DT, logs.SDT, logs.RhoB, rho_fl, k_fl, rho_g, k_g, k0, logs.Porosity)
Vsgc = Vs_GreenbergCastagna(logs.DT,logs.RhoB,k_qz,mu_qz,k_fl,k_b,logs.PhieTotal,rho_fl)
vpggc, vsggc, rhogasgc, kggc = frm(logs.DT, Vsgc, logs.RhoB, rho_fl, k_fl, rho_g, k_g, k0, logs.Porosity)
logs['VS_GC']= Vsgc*1000

sand_cutoff = 0.5
brine_sand = ((logs.Vshale <= sand_cutoff) & (logs.SW >= 0.9))
gas_sand = ((logs.Vshale <= sand_cutoff) & (logs.SW < 0.9))
shale = (logs.Vshale > sand_cutoff)

logs['VP']= (1/(logs.DT*10**-3)*0.3)*10**3
logs['VS'] = (1/(logs.SDT*10**-3)*0.3)*10**3
logs['RHO'] =logs.RhoB

logs['VP_FRMG'] = (1/(logs.DT*10**-3)*0.3)*10**3
logs['VS_FRMG'] = (1/(logs.SDT*10**-3)*0.3)*10**3
logs['RHO_FRMG'] = logs.RhoB
logs['VP_FRMG'][brine_sand|gas_sand] = vpg[brine_sand|gas_sand]
logs['VS_FRMG'][brine_sand|gas_sand] = vsg[brine_sand|gas_sand]
logs['RHO_FRMG'][brine_sand|gas_sand] = rhogas[brine_sand|gas_sand]
logs['IP_FRMG'] = logs.VP_FRMG*logs.RHO_FRMG
logs['IP_FRMG'][shale]= logs.IP
logs['VPVS_FRMG'] = logs.VP_FRMG/logs.VS_FRMG

logs['VP_FRMG_gc'] = (1/(logs.DT*10**-3)*0.3)*10**3
logs['VS_FRMG_gc'] = logs.VS_GC
logs['RHO_FRMG_gc'] = logs.RhoB
logs['VP_FRMG_gc'][brine_sand|gas_sand] = vpg[brine_sand|gas_sand]
logs['RHO_FRMG_gc'][brine_sand|gas_sand] = rhogasgc[brine_sand|gas_sand]
logs['IP_FRMG_gc'] = logs.VP_FRMG_gc*logs.RHO_FRMG_gc
logs['IP_FRMG_gc'][shale]= logs.IP
logs['VPVS_FRMG_gc'] = logs.VP_FRMG_gc/logs.VS_FRMG_gc

Imp1 = logs.IP.values
Imp2= logs.IP_FRMG_gc.values
Imp3 = logs.IP_FRMG.values
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
    f=40            #wavelet frequency
    length=0.512    #Wavelet vector length
    dt=dt           # Sampling prefer to use smiliar to resampled AI
    t0, w = ricker(f, length, dt) # ricker wavelet
    synthetic = np.convolve(w, Rc_tdom, mode='same')
    return synthetic

logs['SS_Original'] = Sintetic_Seismogram(Imp1,t)
logs['SS_Gasgreenbergcastagna'] = Sintetic_Seismogram(Imp2,t)
logs['SS_Gas'] = Sintetic_Seismogram(Imp3,t)

ztop = 6740; zbot = 6800
ll = logs[(logs.Depth >= ztop) & (logs.Depth <= zbot)]


f, ax = plt.subplots(nrows=1, ncols=9, figsize=(10, 80))
f.suptitle('Fluid Replacement Model : WELL-TAUK-1')
ax[0].plot(ll.Vshale, ll.Depth, '-g', label='Vsh')
ax[0].plot(ll.SW, ll.Depth, '-b', label='Sw')
ax[0].plot(ll.PhieTotal,ll.Depth,'-r',label='Fraction Rho')
ax[0].plot(ll.NPHI, ll.Depth, '-k', label='phi')
y2 = ll.Vshale
y1 = y2 * 0 + 0.5

ax[0].fill_betweenx(ll.Depth, y1, y2, where=(y1 >= y2), color='gold', linewidth=0)
ax[0].fill_betweenx(ll.Depth, y1, y2, where=(y1 <= y2), color='lime', linewidth=0)
ax[0].fill_betweenx(ll.Depth, ll.PhieTotal,ll.NPHI , where=(ll.PhieTotal >= ll.NPHI), color='0.5')
ax[0].fill_betweenx(ll.Depth, ll.PhieTotal, ll.NPHI, where=(ll.PhieTotal <= ll.NPHI), color='orange')

ax[1].plot(ll.VP_FRMG, ll.Depth, '-r',label='Gas')
ax[1].plot(ll.VP, ll.Depth, '-', color='0.5',label='Original')


ax[2].plot(ll.VS_FRMG, ll.Depth, '-r',label='Gas')
ax[2].plot(ll.VS_GC, ll.Depth, '-b',label='Original_GC')
ax[2].plot(ll.VS, ll.Depth, '-', color='0.5',label='Original')

ax[3].plot(ll.RHO_FRMG, ll.Depth, '-r',label='Gas')
ax[3].plot(ll.RHO_FRMG_gc, ll.Depth,color='0.6',label='Gas_GC')
ax[3].plot(ll.RHO,ll.Depth, '-',color='0.5',label='Original')

ax[4].plot(ll.IP_FRMG, ll.Depth, '-r',label='Gas')
ax[4].plot(ll.IP_FRMG_gc, ll.Depth, color='0.6',label='Gas_GC')
ax[4].plot(ll.IP, ll.Depth, '-', color='0.5',label='Original')

ax[5].plot(ll.VPVS_FRMG, ll.Depth, '-r',label='Gas')
ax[5].plot(ll.VPVS_FRMG_gc, ll.Depth, color='0.6',label='Gas_GC')
ax[5].plot(ll.Vp0Vs, ll.Depth, '-', color='0.5',label='Original')
for i in range(4):
    ax[6].plot(ll.SS_Original+i,ll.Depth, '-k',label='Original')
    ax[6].fill_betweenx(ll.Depth,ll.SS_Original+i, i,  ll.SS_Original+i > i,  color='k')
    ax[7].plot(ll.SS_Gasgreenbergcastagna+i, ll.Depth, '-b',label='GreenbergCastagna')
    ax[7].fill_betweenx(ll.Depth,ll.SS_Gasgreenbergcastagna+i, i,  ll.SS_Gasgreenbergcastagna+i > i,  color='b')
    ax[8].plot(ll.SS_Gas+i, ll.Depth, '-r',label='Gas')
    ax[8].fill_betweenx(ll.Depth,ll.SS_Gas+i, i,  ll.SS_Gas+i > i,  color='r',)
for i in range(9):
    ax[i].set_ylim(ztop,zbot)
    ax[i].invert_yaxis()
    ax[i].grid()
    ax[i].locator_params(axis='x', nbins=2)
ax[0].legend(fontsize='small', loc='lower right')
ax[0].set_xlabel("Vcl/phi/rho/Sw")
ax[1].legend(fontsize='small', loc='lower right')
ax[1].set_xlabel("Vp [m/s]")
ax[2].legend(fontsize='small', loc='lower right')
ax[2].set_xlabel("Vs[m/s]")
ax[3].legend(fontsize='small', loc='lower right')
ax[3].set_xlabel('Rho [gr/cc]')
ax[4].legend(fontsize='small', loc='lower right')
ax[4].set_xlabel("AI")
ax[5].legend(fontsize='small', loc='lower right')
ax[5].set_xlabel("Vp/Vs")
ax[6].set_xlabel('Original')
ax[7].set_xlabel(" Greenberg Castagna")
ax[8].set_xlabel("gas")
ax[1].set_yticklabels([])
ax[2].set_yticklabels([])
ax[3].set_yticklabels([])
ax[4].set_yticklabels([])
ax[5].set_yticklabels([])
ax[6].set_yticklabels([])
ax[7].set_yticklabels([])
ax[8].set_yticklabels([])

plt.show()