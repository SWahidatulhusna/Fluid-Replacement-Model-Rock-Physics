

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
logs = pd.read_csv('welltaukfacies_lfc.csv')
ll = logs[(logs.Depth >= 6728) & (logs.Depth <= 6728)]
print(ll.PhieTotal)
def msat(vs, rho):
    '''
    Miu saturated

    Input:
    vs = Shear Velocity (Km/s)
    rho = Density (g/cc)

    '''
    return rho * (vs ** 2)


def ksat1(vp, vs, rho):
    '''
    K saturated

    Input:
    vp = P wave velocity (Km/s)
    vs = Shear wave Velocity (Km/s)
    rho = Density (g/cc)

    '''
    Ksat1 = ((vp ** 2) * rho) - 4 / 3 * (msat(vs, rho))
    return np.round(Ksat1, 3)


def kfl(sw, khc, kw):
    '''
    K Fluid

    Input:
    sw = Water Saturation
    khc = K Hydrocarbon
    kw = K Water

    '''
    kfl = ((sw / kw) + ((1 - sw) / khc)) ** -1
    return np.round(kfl, 3)


def km(por, fq, fc, kq, kc):
    '''
    K mineral (K0) 

    Input:
    por = Porosity
    fq = Quartz Fraction
    fc = Calcite Fraction
    kq = Quarzt Bulk
    kc = Calcite Bulk

    '''
    kv = ((fq) * kq) + ((fc) * kc)

    kr = (((fq) / kq) + ((fc) / kc)) ** (-1)

    km = (kv + kr) / 2

    return np.round(km, 3)


def ksat2(ksat1, k0, kfl1, kfl2, por):
    '''
    K Saturation after fluid substitution (C) hpratama 2020

    Input:
    ksat1 = K Saturation before Fluid Substitution
    k0 = K Mineral
    kfl1 = K Fluid before Substitution
    kfl2 = K Fluid after Substitution
    por = Porosity

    '''
    a = (ksat1 / (k0 - ksat1)) - (kfl1 / (por * (k0 - kfl1))) + (kfl2 / (por * (k0 - kfl2)))
    ksat2 = (a * (k0)) / (1 + a)
    return np.round(ksat2, 3)


def rhofluid(sw, so, rhoo, rhow):
    rho1 = ((1 - so) * rhow) + (so * rhoo)
    rho2 = ((sw) * rhow) + ((1 - sw) * rhoo)
    rho1 = np.round(rho1, 3)
    rho2 = np.round(rho2, 3)
    return rho1, rho2


def rhosat(rho, por, rhow, sw, so, rhoo):
    rho1, rho2 = rhofluid(sw, so, rhoo, rhow)

    rhosat = rho + (por * (rho2 - rho1))
    return np.round(rhosat, 3)

def vp(ksat2, vs, rhosat):
    vp = (((ksat2 + (4/3*(msat(vs, rhosat)))) / (rhosat))**0.5)*10**3
    return np.round(vp, 3)
Vp = 3.525
Vs = 1.835

Sw = np.linspace(0,1,51)
So = 1
# So2 = np.linspace(1,0,51)
# Sw2 = 0
Por = 0.093575758
SaturOil= np.linspace(1,0,51)
Fc = 0.3
Fq = 0.7

Kq = 37
Kc = 75
Kg = 0.02
Kw = 2.6

Rho = 2.65
Rhog = 0.1
Rhow = 1.1

# To get new value of Vp, there is some variable that need to calculate.
ksatfl2=[]
# First we need to calculate K saturation before fluid subtitution.
# This K saturartion could be calculated from the old
Ksat1 = ksat1(Vp, Vs, Rho)
# Then, we estimate the K fluid before and after fluid substitution
Kfl1 = kfl((1-So), Kg, Kw)
for i in range (len(Sw)):
    Kfl2 = kfl(Sw[i], Kg, Kw)
    ksatfl2.append(Kfl2)
# After that, we evaluate K Mineral that can be
# estimate from Hill formula, this formula need Voight-Reuss Bounds calculation too.
# First we need to calculate K saturation before fluid subtitution.
K0 = km(Por, Fq, Fc, Kq, Kc)
# print(ksatfl2)
# Here we calculate K saturated after fluid substitution
ksatur2=[]
for j in range (len(ksatfl2)):
    Ksat2 = ksat2(Ksat1, K0, Kfl1, ksatfl2[j], Por)
    ksatur2.append(Ksat2)
print(ksatur2)
# To calculate the new Vp, we need to calculate density saturation.
# After that we estimate new Vp with new K saturated after fluid substitution
rhosatur=[]
for k in range (len(ksatur2)):
    Rhosat = rhosat(Rho, Por, Rhow, Sw[k], So, Rhog)
    rhosatur.append(Rhosat)
# print(rhosatur)

Velocityprimer=[]
for l in range (len(rhosatur)):
    Vp = vp(ksatur2[l], Vs, rhosatur[l])
    Velocityprimer.append(Vp)

# ksatfl21=[]
# # print(Velocityprimer)
# Kfl21 = kfl((1-Sw2), Kg, Kw)
# for i in range (len(So2)):
#     Kfl21 = kfl(So2[i], Kg, Kw)
#     ksatfl21.append(Kfl21)
# # After that, we evaluate K Mineral that can be
# # estimate from Hill formula, this formula need Voight-Reuss Bounds calculation too.
# # First we need to calculate K saturation before fluid subtitution.
# K0 = km(Por, Fq, Fc, Kq, Kc)
# # print(ksatfl2)
# # Here we calculate K saturated after fluid substitution
# ksatur21=[]
# for j in range (len(ksatfl2)):
#     Ksat21 = ksat2(Ksat1, K0, Kfl1, ksatfl21[j], Por)
#     ksatur21.append(Ksat21)
# # print(ksatur2)
# # To calculate the new Vp, we need to calculate density saturation.
# # After that we estimate new Vp with new K saturated after fluid substitution
# rhosatur1=[]
# for k in range (len(ksatur2)):
#     Rhosat = rhosat(Rho, Por, Rhow, Sw[k], So, Rhog)
#     rhosatur1.append(Rhosat)
# # print(rhosatur)
#
# Velocityprimer2=[]
# for l in range (len(rhosatur1)):
#     Vp = vp(ksatur2[l], Vs, rhosatur1[l])
#     Velocityprimer2.append(Vp)




fig,ax = plt.subplots()


ax.plot(Sw,Velocityprimer,'k-')
# ax.plot(SaturOil,Velocityprimer2  ,'k-')

plt.show()