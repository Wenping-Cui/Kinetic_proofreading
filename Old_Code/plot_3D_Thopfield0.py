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
fig_name='3D_Thopfield3';
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

KT = 1.0*10**(4);    # Kinase total
sigma = 0;   # Kinase phosphorylation rate
epsilon = 0.01; # Kinase dephosphorylation rate
K = KT/(1.0 + sigma/epsilon*(C[m]+D[m]));  # Kinase


alpha = 3.0*10**(-4); # alpha*K  the last forward rate 
b = 0.1;  # Back rate
phi = 1.0; # Forward rate
gamma = 10**(-3); 
W = 100.0; # Reaction rate to the final product
# Summary the parameters
Num = 10**5; # Number of scatter points
def data_3Dhopfied(Num):
    data = np.zeros((Num,8));
    for num in range(0, Num):    
        phi = 10**np.random.uniform(-10,8);
        b = 10**np.random.uniform(-12,12);
        gamma = 10**np.random.uniform(-4,1);
        alpha = phi/KT
        W=100.0
        para = np.array([Ls,Lf,taus,tauf,R,kappa,KT,sigma,epsilon,alpha,b,phi,gamma,W])
        output = model(N,m,para);
        data[num] = np.concatenate((output[0:3],[output[4][N],output[5][N],b,phi,gamma]), axis=0)
    return data[:,0:8]
data = data_3Dhopfied(Num)
data =data[np.where( data[:,0] > 0.0 )]   
data =data[np.where( data[:,1] > 0.0 )] 
data =data[np.where( data[:,2] < 100.0 )] 
C = data[:,3]
D = data[:,4]
A = C+D
Speed = data[:,0]
Energy = data[:,1]
#Energy = np.log(data[:,1])
Error = data[:,2];


color_bar='rainbow'




fig, ax = plt.subplots()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 23})
my_cmap = matplotlib.cm.get_cmap(color_bar)
sc = plt.scatter(Speed, Energy, c=C,edgecolors='none', vmin=10**(-12), vmax=10**(4), cmap=my_cmap, s=10,norm=matplotlib.colors.LogNorm())
clb = plt.colorbar(sc)
clb.set_label('C', fontsize=15, labelpad=-40, y=-0.028, rotation=0)
plt.xlabel('MFPT$(sec)$', fontsize=25)
plt.ylabel('Dissipation$(k_B T)$', fontsize=25)
ax.set_xlim([10**(0),10**13])
ax.set_ylim([10**(1),10**6])
ax.set_xscale("log")
ax.set_yscale("log")
fig_name='3D_Thopfield_C';
fig.set_size_inches(14, 12)
fig.savefig(fig_name, transparent=True, bbox_inches='tight', \
                        pad_inches=0)




fig, ax = plt.subplots()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 23})
my_cmap = matplotlib.cm.get_cmap(color_bar)
sc = plt.scatter(Speed, Energy, c=D,edgecolors='none', vmin=10**(-12), vmax=10**(4), cmap=my_cmap, s=10,norm=matplotlib.colors.LogNorm())
clb = plt.colorbar(sc)
clb.set_label('D', fontsize=15, labelpad=-40, y=-0.028, rotation=0)
plt.xlabel('MFPT$(sec)$', fontsize=25)
plt.ylabel('Dissipation$(k_B T)$', fontsize=25)
ax.set_xlim([10**(0),10**13])
ax.set_ylim([10**(1),10**6])
ax.set_xscale("log")
ax.set_yscale("log")
fig_name='3D_Thopfield_D';
fig.set_size_inches(14, 12)
fig.savefig(fig_name, transparent=True, bbox_inches='tight', \
                        pad_inches=0)

fig, ax = plt.subplots()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 23})
my_cmap = matplotlib.cm.get_cmap(color_bar)
sd = plt.scatter(Speed, Energy, c=A,edgecolors='none', vmin=10**(-12), vmax=10**(4), cmap=my_cmap, s=10,norm=matplotlib.colors.LogNorm())
clb = plt.colorbar(sd)
clb.set_label('Output Signal', fontsize=15, labelpad=-40, y=-0.028, rotation=0)
plt.xlabel('Time$(sec)$', fontsize=25)
plt.ylabel('Dissipation($k_B T)$', fontsize=25)
ax.set_xlim([10**(0),10**13])
ax.set_ylim([10**(0),10**6])
ax.set_xscale("log")
ax.set_yscale("log")
fig_name='3D_Thopfield_Output';
fig.set_size_inches(14, 12)
fig.savefig(fig_name, transparent=True, bbox_inches='tight', \
                        pad_inches=0)

fig, ax = plt.subplots()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 23})
my_cmap = matplotlib.cm.get_cmap(color_bar)
se = plt.scatter(Speed, Energy, c=Error,edgecolors='none', vmin=10**(-4), vmax=10**(2), cmap=my_cmap, s=10,norm=matplotlib.colors.LogNorm())
clb = plt.colorbar(se)
clb.set_label('Error rate', fontsize=15, labelpad=-40, y=-0.028, rotation=0)
plt.xlabel('MFPT$(sec)$', fontsize=25)
plt.ylabel('Dissipation($k_B T)$', fontsize=25)
ax.set_xlim([10**(0),10**13])
ax.set_ylim([10**(0),10**6])
ax.set_xscale("log")
ax.set_yscale("log")
fig_name='3D_Thopfield';
fig.set_size_inches(14, 12)
fig.savefig(fig_name, transparent=True, bbox_inches='tight', \
                        pad_inches=0)

plt.subplot(1, 2, 1)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 17})
my_cmap = matplotlib.cm.get_cmap(color_bar)
sc = plt.scatter(Speed, Energy, c=Error,edgecolors='none', cmap=my_cmap, s=10,norm=matplotlib.colors.LogNorm())
#sc = plt.scatter(Speed, Energy, c=Error,edgecolors='none', vmin=10**(-4), vmax=10**(4), cmap=my_cmap, s=10,norm=matplotlib.colors.LogNorm())
clb = plt.colorbar(sc,orientation="horizontal")
clb.set_label('Error rate', fontsize=18, labelpad=-36, y=1.2, rotation=0)
#plt.text(-0.05, 1.1, 'a', transform=plt.gca().transAxes,
#      fontsize=25, fontweight='bold', va='top', ha='right')
plt.xlabel('MFPT$(sec)$', fontsize=25)
plt.ylabel('Dissipation($k_B T/s)$', fontsize=25)
#plt.axvline(10**2, color='k', linestyle='--',linewidth=3.0)
plt.xlim([10**(0),10**13])
plt.ylim([10**(1),10**6])
plt.xscale("log")
plt.yscale("log")
plt.subplot(1, 2, 2)
my_cmap = matplotlib.cm.get_cmap(color_bar)
sd = plt.scatter(Speed, Energy, c=A,edgecolors='none', vmin=10**(-12), vmax=10**(4), cmap=my_cmap, s=10,norm=matplotlib.colors.LogNorm())
clb = plt.colorbar(sd,orientation="horizontal")
clb.set_label('Output Signal', fontsize=18, labelpad=-36, y=1.2, rotation=0)
plt.xlabel('MFPT$(sec)$', fontsize=25)
plt.ylabel('Dissipation($k_B T/s)$', fontsize=25)
plt.xlim([10**(0),10**13])
plt.ylim([10**(1),10**6])
#plt.text(-0.05, 1.1, 'b', transform=plt.gca().transAxes,
#      fontsize=25, fontweight='bold', va='top', ha='right')
#plt.axvline(10**2, color='k', linestyle='--',linewidth=3.0)
plt.xscale("log")
plt.yscale("log")
fig_name='3D_Thopfield_2';
fig.set_size_inches(14, 7)
plt.savefig(fig_name, transparent=True, bbox_inches='tight', \
                        pad_inches=0)




