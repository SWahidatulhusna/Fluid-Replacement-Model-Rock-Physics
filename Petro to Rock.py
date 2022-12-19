import matplotlib.pyplot as plt
import lasio
import pandas as pd
import numpy as np
las = lasio.read("WELL_TAUK-1.las")
df = las.df()
df.describe()
rho_ma = 2.64
rho_f = 1.08
dt_ma = 55.5
dt_f = 189

file = './WELL_TAUK-1.las'
data1 = np.loadtxt(file,skiprows=38)


memories = ['Depth','Checkshot','RhoB','GR','NPHI','DT','PHIE_1','PHIT_1','SDT','SW']
data1 = pd.DataFrame(data1,columns=memories)
data1 = data1[['Depth','Checkshot','RhoB','GR','NPHI','DT','PHIE_1','PHIT_1','SDT','SW']]
sw = np.linspace(0,1,51)
Simulation_SW = [[0 for j in range(0,801)]for i in range(len(sw))]
# print(ST)
for j in range (len(sw)):
    for i in range(0,801):
        Simulation_SW[j][i] = i*0 + sw[j]



# Homework1
V_shale = (las['GR'] - min(las['GR'])) / (max(las['GR']) - min(las['GR']))
Porosity_Total = (rho_ma - las['RHOB']) / (rho_ma - rho_f)
PHI_DT = (las['DT'] - dt_ma) / (dt_f - dt_ma)
Vp = (1 / (las['DT'] * 10 ** -3) * 0.3) * 10 ** 3
Vs = (1 / (las['SDT'] * 10 ** -3) * 0.3) * 10 ** 3
PHIE = Porosity_Total * (1 - V_shale)
Porosity = (PHIE - min(PHIE)) / (max(PHIE) - min(PHIE))
ps = Vp / Vs
IP= Vp * las['RHOB']
data1['Vshale'] = V_shale
data1['PhieTotal'] = Porosity_Total
data1['PhieDT'] = PHI_DT
data1['PhiE'] = PHIE
data1['Porosity'] = Porosity
data1['Vp0Vs'] = ps
data1['SW100'] = Simulation_SW[50]
data1['SW98'] = Simulation_SW[49]
data1['SW96'] = Simulation_SW[48]
data1['SW94'] = Simulation_SW[47]
data1['SW92'] = Simulation_SW[46]
data1['SW90'] = Simulation_SW[45]
data1['SW88'] = Simulation_SW[44]
data1['SW86'] = Simulation_SW[43]
data1['SW84'] = Simulation_SW[42]
data1['SW82'] = Simulation_SW[41]
data1['SW80'] = Simulation_SW[40]
data1['SW78'] = Simulation_SW[39]
data1['SW76'] = Simulation_SW[38]
data1['SW74'] = Simulation_SW[37]
data1['SW72'] = Simulation_SW[36]
data1['SW70'] = Simulation_SW[35]
data1['SW68'] = Simulation_SW[34]
data1['SW66'] = Simulation_SW[33]
data1['SW64'] = Simulation_SW[32]
data1['SW62'] = Simulation_SW[31]
data1['SW60'] = Simulation_SW[30]
data1['SW58'] = Simulation_SW[29]
data1['SW56'] = Simulation_SW[28]
data1['SW54'] = Simulation_SW[27]
data1['SW52'] = Simulation_SW[26]
data1['SW50'] = Simulation_SW[25]
data1['SW48'] = Simulation_SW[24]
data1['SW46'] = Simulation_SW[23]
data1['SW44'] = Simulation_SW[22]
data1['SW42'] = Simulation_SW[21]
data1['SW40'] = Simulation_SW[20]
data1['SW38'] = Simulation_SW[19]
data1['SW36'] = Simulation_SW[18]
data1['SW34'] = Simulation_SW[17]
data1['SW32'] = Simulation_SW[16]
data1['SW30'] = Simulation_SW[15]
data1['SW28'] = Simulation_SW[14]
data1['SW26'] = Simulation_SW[13]
data1['SW24'] = Simulation_SW[12]
data1['SW22'] = Simulation_SW[11]
data1['SW20'] = Simulation_SW[10]
data1['SW18'] = Simulation_SW[9]
data1['SW16'] = Simulation_SW[8]
data1['SW14'] = Simulation_SW[7]
data1['SW12'] = Simulation_SW[6]
data1['SW10'] = Simulation_SW[5]
data1['SW8'] = Simulation_SW[4]
data1['SW6'] = Simulation_SW[3]
data1['SW4'] = Simulation_SW[2]
data1['SW2'] = Simulation_SW[1]
data1['SW0'] = Simulation_SW[0]
data1['IP']= IP



data1.to_csv('Well_Tauk_AfterProcessing.csv',index=False)
# Assumtion
K1 = 37
K2 = 2.6
M1 = 44
M2 = 0
Vp1 = 1.5
Vp2 = 6.5

# Therotical Bound
K_Voigt = Porosity * K2 + (1 - Porosity) * K1
K_Reuss = 1 / (Porosity / K2 + (1 - Porosity) / K1)
K_VRH = (K_Voigt + K_Reuss) / 2
M_Voigt = Porosity * M2 + (1 - Porosity) * M1
M_Reuss = 1 / (Porosity / M2 + (1 - Porosity) / M1)
M_VRH = (M_Voigt + M_Reuss) / 2
Kupper_HS = K1 + Porosity / ((K2 - K1) ** -1 + (1 - Porosity) * (K1 + 1.25 * M1) ** -1)
Klower_HS = K1 + Porosity / ((K2 - K1) ** -1 + (1 - Porosity) * (K1 + 1.25 * M2) ** -1)
Mupper_HS = M1 + Porosity / (((M2 - M1) ** -1) + (2 * (1 - Porosity)) * (K2 + 2 * M1) / (5 * M1 * (K2 + 1.25 * M1)))
Mlower_HS = M1 + Porosity / (((M2 - M1) ** -1) + (2 * (1 - Porosity)) * (K1 + 2 * M1) / (5 * M1 * (K2 + 1.25 * M1)))
V_Voigt = Porosity * Vp1 + (1 - Porosity) * Vp2
V_Reuss = 1 / (Porosity / Vp1 + (1 - Porosity) / Vp2)

# Calculated k,mu,
mcal = (Vs ** 2) * Porosity_Total * 10 ** -5
kcal = ((Vp ** 2) * Porosity_Total - (4 / 3) * mcal) * 10 ** -5


# Homework2
rows, cols = 1, 5
fig, ax = plt.subplots(nrows=rows, ncols=cols, figsize=(12, 10), sharey=True)
# Plot Vshale vs depth
ax[0].plot(V_shale, las['DEPTH'])
ax[0].set_ylim(max(las['DEPTH']), min(las['DEPTH']))
ax[0].minorticks_on()
ax[0].grid(which='major', linestyle='-', linewidth='0.5', color='green')
ax[0].grid(which='minor', linestyle=':', linewidth='0.5', color='black')
ax[0].set_title('Vshale')
y2 = V_shale
y1 = y2 * 0 + 0.5
ax[0].fill_betweenx(las['DEPTH'], y1, y2, where=(y1 >= y2), color='gold', linewidth=0)
ax[0].fill_betweenx(las['DEPTH'], y1, y2, where=(y1 <= y2), color='lime', linewidth=0)
# Plot Porosity total vs depth
ax[1].plot(Porosity_Total, las['DEPTH'])
ax[1].plot(las['NPHI'], las['DEPTH'])
ax[1].set_ylim(max(las['DEPTH']), min(las['DEPTH']))
ax[1].minorticks_on()
ax[1].grid(which='major', linestyle='-', linewidth='0.5', color='green')
ax[1].grid(which='minor', linestyle=':', linewidth='0.5', color='black')
ax[1].set_title('Porositas Total & NPhi')
ax[1].fill_betweenx(las['DEPTH'], Porosity_Total, las['NPHI'], where=(Porosity_Total >= las['NPHI']), color='red')
ax[1].fill_betweenx(las['DEPTH'], Porosity_Total, las['NPHI'], where=(Porosity_Total <= las['NPHI']), color='blue')
# Plot Porosity total sonic vs depth
ax[2].plot(PHI_DT, las['DEPTH'])
ax[2].set_ylim(max(las['DEPTH']), min(las['DEPTH']))
ax[2].minorticks_on()
ax[2].grid(which='major', linestyle='-', linewidth='0.5', color='green')
ax[2].grid(which='minor', linestyle=':', linewidth='0.5', color='black')
ax[2].set_title('PHI_DT')
# Plot Porosity effective vs depth
ax[3].plot(PHIE, las['DEPTH'])
ax[3].set_ylim(max(las['DEPTH']), min(las['DEPTH']))
ax[3].minorticks_on()
ax[3].grid(which='major', linestyle='-', linewidth='0.5', color='green')
ax[3].grid(which='minor', linestyle=':', linewidth='0.5', color='black')
ax[3].set_title('PhiE')
# Plot S Wave vs depth
ax[4].plot(las['SDT'], las['DEPTH'])
ax[4].set_ylim(max(las['DEPTH']), min(las['DEPTH']))
ax[4].minorticks_on()
ax[4].grid(which='major', linestyle='-', linewidth='0.5', color='green')
ax[4].grid(which='minor', linestyle=':', linewidth='0.5', color='black')
ax[4].set_title('SDT')
plt.subplots_adjust(wspace=0.015)

# Homework2
ig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))

# plot Vp/Vs vs porosity
ax[0, 0].scatter(ps, Porosity, c=V_shale, cmap='viridis', label='$Vp/Vs$')
ax[0, 0].grid(which='major', linestyle='-', linewidth='0.2', color='green')
ax[0, 0].set_xlabel("$Vp/vs$")
ax[0, 0].set_ylabel("$Porosity$")
ax[0, 0].set_facecolor('#90EE90')

# plot K vs porosity
ax[1, 0].plot(Porosity, K_Voigt, linewidth=2, color='blue', label='$Voight$')
ax[1, 0].plot(Porosity, K_Reuss, linewidth=2, color='red', label='$Reuss$')
ax[1, 0].plot(Porosity, K_VRH, linewidth=2, color='purple', label='$Hill$')
ax[1, 0].plot(Porosity, Kupper_HS, linewidth=2, color='yellow', label='$Hashin Strihksman+$')
ax[1, 0].plot(Porosity, Klower_HS, linewidth=1, color='orange', label='$Hashin Strihksman-$')
ax[1, 0].scatter(Porosity, kcal, c=V_shale, cmap='viridis')
ax[1, 0].legend(loc='upper center', bbox_to_anchor=(0.3, 1.05), prop={'size': 5})
ax[1, 0].set_xlabel("$Porosity$")
ax[1, 0].set_ylabel("$Bulk Modulus(Gpa)$")
ax[1, 0].grid(which='major', linestyle='-', linewidth='0.2', color='green')
ax[1, 0].set_facecolor('#90EE90')

# plot M vs porosity
ax[0, 1].plot(Porosity, M_Voigt, linewidth=2, color='blue', label='$Voight$')
ax[0, 1].plot(Porosity, M_Reuss, linewidth=2, color='red', label='$Reuss$')
ax[0, 1].plot(Porosity, M_VRH, linewidth=2, color='purple', label='$Hill$')
ax[0, 1].plot(Porosity, Mupper_HS, linewidth=2, color='yellow', label='$Hashin Strihksman+$')
ax[0, 1].plot(Porosity, Mlower_HS, linewidth=2, color='orange', label='$Hashin Strihksman-$')
ax[0, 1].scatter(Porosity, mcal, c=V_shale, cmap='viridis')
ax[0, 1].legend(bbox_to_anchor=(0.3, 1.05),prop={'size': 5})
ax[0, 1].set_xlabel("$Porosity$")
ax[0, 1].set_ylabel("$Shear Modulus(Gpa)$")
ax[0, 1].grid(which='major', linestyle='-', linewidth='0.2', color='green')
ax[0, 1].set_facecolor('#90EE90')

# plot Vp vs porosity
ax[1, 1].plot(Porosity, V_Voigt * 10 ** 3, linewidth=2, color='blue', label='$Voight$')
ax[1, 1].plot(Porosity, V_Reuss * 10 ** 3, linewidth=2, color='red', label='$Reuss$')
ax[1, 1].scatter(Porosity, Vp, c=V_shale, cmap='viridis')
ax[1, 1].grid(which='major', linestyle='-', linewidth='0.2', color='green')
ax[1, 1].legend(bbox_to_anchor=(1.1, 1.05))
ax[1, 1].set_xlabel("$Porosity$")
ax[1, 1].set_ylabel("$Vp (m/s)$")
ax[1, 1].set_facecolor('#90EE90')

plt.subplots_adjust(hspace=0.23)
plt.subplots_adjust(hspace=0.3)
cax = ig.add_axes([0.925, 0.1, 0.02, 0.8])
ig.colorbar(ax[0, 1].scatter(Porosity, mcal, c=V_shale, cmap='viridis'), cax=cax)
plt.show()
