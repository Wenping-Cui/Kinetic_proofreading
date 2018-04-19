import time
import os
import math
import numpy as np
from Adaption_lib1 import *
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import colors
# These are the "Tableau 20" colors as RGB.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)
fig_name='2D_Thopfield3';
N = 4;
m = 2;
# Definition of parameters
C = np.zeros(N+1);  # Agonist complex phosphorylated n times
D = np.zeros(N+1);  # Non-Agonist complex phosphorylated n times

Ls = 10**4;    # Concentration of self ligands 10^4 - 10^5
Lf = 10**4;       # Concentration of foreign ligands ~ 10 per cell
taus = 1.0;      # Self complex bonding time
tauf = 10.0;     # Foreign complex bonding time


R = 10.0**4;     # Concentration of receptor when Ls = 0
kappa = 3.0*10**(2); # Receptor-ligand bonding rate
kappaR = 3.0;    # Considering Ls, suppose the receptors are largely in excess

KT = 1.0*10**(3);    # Kinase total 
sigma = 0;   # Kinase phosphorylation rate
epsilon = 0.01; # Kinase dephosphorylation rate
K = KT/(1.0 + sigma/epsilon*(C[m]+D[m]));  # Kinase


alpha = 3.0*10**(-4); # alpha*K  the last forward rate 
b = 0.1;  # Back rate
phi = 1.0; # Forward rate
gamma = 10**(-3);
W = 100.0; # Reaction rate to the final product
# Summary the parameters

Num = 80000
def data_3Dhopfied(Num):
    data = np.zeros((Num,9));
    for num in range(0, Num):
#        gamma = 10**np.random.uniform(-5,1);
        gamma = 10**(-3)
        phi = 10**np.random.uniform(-10,8);
        b = 10**np.random.uniform(-12,10);
        alpha = phi/KT
        W=100.0
        para = np.array([Ls,Lf,taus,tauf,R,kappa,KT,sigma,epsilon,alpha,b,phi,gamma,W])
        output = model(N,m,para);
        data[num] = np.concatenate((output[0:3],[output[4][N],output[5][N],phi,phi/b,b,gamma]), axis=0)
    return data[:,0:9]
data = data_3Dhopfied(Num)

Num = 10**5
def data_3Dhopfied0(Num):
    data0 = np.zeros((Num,9));
    for num in range(0, Num):    
        gamma = 10**np.random.uniform(-4,1);
        phi = 10**np.random.uniform(-10,8);
        b = 10**np.random.uniform(-12,10);
        alpha = phi/KT
        W=100.0
        para = np.array([Ls,Lf,taus,tauf,R,kappa,KT,sigma,epsilon,alpha,b,phi,gamma,W])
        output = model(N,m,para);
        data0[num] = np.concatenate((output[0:3],[output[4][N],output[5][N],phi,phi/b,b,gamma]), axis=0)
    return data0[:,0:9]
data0 = data_3Dhopfied0(Num)





data2 =data0[np.where( data0[:,0] > 10**6.5)]
data2 =data2[np.where( data2[:,0] < 10.0**7.5)]
data2 =data2[np.where( data2[:,1] > 0.0 )]
#data2 =data2[np.where( data2[:,2] < 1.0)]

Speed2 = data2[:,0]

Energy2 = data2[:,1]

print max(Energy2)
#Energy = np.log(data[:,1])
Error2 = data2[:,2];
Phi2 = data2[:,5];
Phi_b2 = 1/data2[:,6];
B2 = data2[:,7];
Gamma2=data2[:,8]
data3 =data0[np.where( data0[:,1] > 10.0**3.3)]
data3 =data3[np.where( data3[:,1] < 10.0**3.7)]
data3 =data3[np.where( data3[:,0] > 0.0 )]
#data3 =data3[np.where( data3[:,2] < 1.0)]

Speed3 = data3[:,0]
Energy3 = data3[:,1]
#Energy = np.log(data[:,1])
Error3 = data3[:,2];
Phi3 = data3[:,5];
Phi_b3 = 1/data3[:,6];
B3 = data3[:,7];
Gamma3=data3[:,8]


data_2 =data[np.where( data[:,0] > 10**6.9)]
data_2 =data_2[np.where( data_2[:,0] < 10.0**7.1)]
data_2 =data_2[np.where( data_2[:,1] > 0.0 )]
#data2 =data2[np.where( data2[:,2] < 1.0)]

Speed_2 = data_2[:,0]
Energy_2 = data_2[:,1]
#Energy = np.log(data[:,1])
Error_2 = data_2[:,2];
Phi_2 = data_2[:,5];
Phi__b2 = 1/data_2[:,6];
B_2 = data_2[:,7];
Gamma_2=data_2[:,8]
data_3 =data[np.where( data[:,1] > 10.0**3.1)]
data_3 =data_3[np.where( data_3[:,1] < 10.0**3.4)]
data_3 =data_3[np.where( data_3[:,0] > 0.0 )]
#data3 =data3[np.where( data3[:,2] < 1.0)]

Speed_3 = data_3[:,0]
Energy_3 = data_3[:,1]
#Energy = np.log(data[:,1])
Error_3 = data_3[:,2];
Phi3 = data3[:,5];
Phi__b3 = 1/data_3[:,6];
B_3 = data_3[:,7];
Gamma_3=data_3[:,8]

fig, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 23})
ax1.scatter(Energy2, Error2, s=25, facecolors='none', edgecolors=tableau20[4],label=' Random $\\gamma$')
ax1.scatter(Energy_2, Error_2, s=20, facecolors='none', edgecolors=tableau20[6],label=' $\\gamma = 10^{-3}$')
legend = ax1.legend(loc='lower left',frameon=False, fontsize=20)
ax1.set_xlabel('Dissipation$(k_B T/s)$', fontsize=25)
ax1.set_ylabel('Error rate', fontsize=25)
ax1.set_ylim([10**(-4),10**0.5])
ax1.set_xlim([10**(1),10**3.5])
ax1.set_xscale("log")
ax1.set_yscale("log")

ax2.scatter(Speed3, Error3, s=25, facecolors='none', edgecolors=tableau20[4],label=' Random $\\gamma$')
ax2.scatter(Speed_3, Error_3, s=20, facecolors='none', edgecolors=tableau20[6],label=' $\\gamma = 10^{-3}$')
ax2.set_xlabel('MFPT$(sec)$', fontsize=25)
ax2.set_ylabel('Error rate', fontsize=25)
ax2.set_ylim([10**(-4),10**0.5])
ax2.set_xlim([10**(3),10**8.5])
ax2.set_xscale("log")
ax2.set_yscale("log")
legend = ax2.legend(loc='lower right',frameon=False, fontsize=20)
fig_name='2D_Thopfield3';
fig.set_size_inches(18, 7)
fig.savefig(fig_name, transparent=True, bbox_inches='tight',pad_inches=0)
plt.show()
