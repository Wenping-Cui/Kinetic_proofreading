# Code to solve parallel adaptive sorting 
# Specific case: N=4 and m=2
import numpy as np
import matplotlib.pyplot as plt
from Adaption_lib import *
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
N = 4;
m = 2;
# Definition of parameters
C = np.zeros(N+1);  # Agonist complex phosphorylated n times
D = np.zeros(N+1);  # Non-Agonist complex phosphorylated n times

L = 2.0*10**4
Ls = 10**4;    # Concentration of self ligands 10^4 - 10^5
Lf = 10**4;       # Concentration of foreign ligands ~ 10 per cell
taus = 1.0;      # Self complex bonding time
tauf = 10.0;     # Foreign complex bonding time

R = 10.0**4;     # Concentration of receptor when Ls = 0
kappa = 3.0*10**(2); # Receptor-ligand bonding rate
kappaR = 3.0;    # Considering Ls, suppose the receptors are largely in excess

KT = 1.0*10**(4);    # Kinase total
KT = 1.0*10**(4);    # Kinase total
sigma = 1.0;   # Kinase phosphorylation rate
epsilon = 2.0; # Kinase dephosphorylation rate
K = KT/(1.0 + sigma/epsilon*(C[m]+D[m]));  # Kinase


alpha = 3.0*10**(-4); # alpha*K  the last forward rate
b = 0.1;  # Back rate
phi = 1.0; # Forward rate
gamma = 10**(-3); 
W = 100.0; # Reaction rate to the final product
# Summary the parameters

para = np.array([Ls,Lf,taus,tauf,R,kappa,KT,sigma,epsilon,alpha,b,phi,gamma,W])

# Output: speed, energy and eta
def data_ana(p,sigma):
    Lf = p*L
    Ls = L-p*L
    Num = 4*10**2; # Number of scatter points
    data = np.zeros((Num,4+6));
    for num in range(0, Num):
#        kappa = np.random.uniform(1.0/tauf/R,10**(2));
#        sigma = 10**np.random.uniform(-2,0.5);
#        epsilon = 10**np.random.uniform(-2,0.5);
#        alpha = 10**np.random.uniform(-3,2); 
        phi = 10**(-4.5+(4.5+1.2)/Num*num);   
#        b = 10**np.random.uniform(-4,-1)*phi;
        b = 10**(-2)*phi;
#        kappa = phi/R;
#        sigma = 0;
#        alpha = phi/KT
        para = np.array([Ls,Lf,taus,tauf,R,kappa,KT,sigma,epsilon,alpha,b,phi,gamma,W])
        output = model(4,2,para);
        data[num] = np.concatenate((output[0:4],[kappa,sigma,epsilon,alpha,b,phi]), axis=0)
    return data[:,0:4+6]
def data_hop(p,sigma):
    Lf = p*L
    Ls = L-p*L
    Num = 4*10**2; # Number of scatter points
    data = np.zeros((Num,4+6));
    sigma=0
    for num in range(0, Num):
#        kappa = 10**np.random.uniform(-4,4);    
#        sigma = 10**np.random.uniform(-2,0.5);
#        epsilon = 10**np.random.uniform(-2,0.5);
#        alpha = 10**np.random.uniform(-4,4);    
        phi = 10**(-4.5+(4.5+1.2)/Num*num);
        b = 0.01*phi;
#        kappa = phi/R;
#        sigma = 0;
        alpha = phi/KT
        para = np.array([Ls,Lf,taus,tauf,R,kappa,KT,sigma,epsilon,alpha,b,phi,gamma,W])
        output = model(N,m,para);
        data[num] = np.concatenate((output[0:4],[kappa,sigma,epsilon,alpha,b,phi]), axis=0)
    return data[:,0:4+6]
# Sorting by eta(error rate)

# Define bin size(of the array size) for energy consumption
bin_error = (taus/tauf)**(N) ;
print bin_error
data1 = data_ana(0.01,sigma);
data2 = data_ana(0.1,sigma);
data3 = data_ana(0.5,sigma);
data4 = data_hop(0.01,sigma);
data5 = data_hop(0.1,sigma);
data6 = data_hop(0.5,sigma);
#data4 = data_ana(0.001);
out = data3[data3[:,2].argsort()]
print out[0][9]
out = np.transpose(out);

data1 = np.transpose(data1);
data2 = np.transpose(data2);
data3 = np.transpose(data3);
data4 = np.transpose(data4);
data5 = np.transpose(data5);
data6 = np.transpose(data6);

f, ((ax1, ax2, ax3)) = plt.subplots(1, 3, sharex=False, sharey=False)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 25})
ax1.text(-0.11, 1.07, 'a', transform=ax1.transAxes,
      fontsize=25, fontweight='bold', va='top', ha='right')
ax2.text(-0.11, 1.07, 'b', transform=ax2.transAxes,
      fontsize=25, fontweight='bold', va='top', ha='right')
ax3.text(-0.11, 1.07, 'c', transform=ax3.transAxes,
      fontsize=25, fontweight='bold', va='top', ha='right')
ax1.plot(data1[0], data1[2],color=tableau20[2],linewidth=3.0, label=' $p = 0.01$')
ax1.plot(data2[0], data2[2],color=tableau20[4],linewidth=3.0,label=' $p = 0.1$')
ax1.plot(data3[0], data3[2],color=tableau20[8],linewidth=3.0,label=' $p = 0.5$')
ax1.plot(data4[0], data4[2], 'r--',color=tableau20[3],linewidth=3.0)
ax1.plot(data5[0], data5[2], 'r--',color=tableau20[5],linewidth=3.0)
ax1.plot(data6[0], data6[2], 'r--',color=tableau20[9],linewidth=3.0,label='KPR')
ax1.set_xscale('log')
ax1.set_yscale('log')
#ax1.scatter(speed_boundary, eta_boundary, s=20, c='r')
ax1.set_ylim([10**(-4),10**3])
ax1.set_xlim([10**(1),10**11])

#ax1.set_xlabel('Time', fontsize=25)
legend = ax1.legend(loc='upper right',frameon=False, fontsize=25)
ax1.set_xlabel('MFPT$(sec)$', fontsize=25)
ax1.set_ylabel('Error rate', fontsize=25)


ax1.axvline(10**2, color='k', linestyle='--',linewidth=3.0)
ax2.axvline(10**2, color='k', linestyle='--',linewidth=3.0)
ax2.set_xlim([10**(1),10**11])
ax2.set_ylim([10**(1),10**6])

ax2.plot(data1[0], data1[1],color=tableau20[2],linewidth=3.0, label=' $L_f = 0.01L_t$')
ax2.plot(data2[0], data2[1],color=tableau20[4],linewidth=3.0,label=' $L_f = 0.1L_t$')
ax2.plot(data3[0], data3[1],color=tableau20[8],linewidth=3.0,label=' $L_f = 0.5L_t$')
ax2.plot(data4[0], data4[1], 'r--',color=tableau20[3],linewidth=3.0)
ax2.plot(data5[0], data5[1], 'r--',color=tableau20[5],linewidth=3.0)
ax2.plot(data6[0], data6[1], 'r--',color=tableau20[9],linewidth=3.0,label='KPR')
#ax2.scatter(out[0][0:n_out], out[1][0:n_out],color = 'r', marker = 'h', s = 40,label='Hopfield')
#ax2.set_xlim([10**3.5,10**8.5])
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('MFPT$(sec)$', fontsize=25)
ax2.set_ylabel('Dissipation$(k_B T/s)$', fontsize=25)

ax3.plot(data1[1], data1[2],color=tableau20[2],linewidth=3.0, label=' $L_f = 0.01L_T$')
ax3.plot(data2[1], data2[2],color=tableau20[4],linewidth=3.0,label=' $L_f = 0.1L_t$')
ax3.plot(data3[1], data3[2],color=tableau20[8],linewidth=3.0,label=' $L_f = 0.5L_t$')
ax3.plot(data4[1], data4[2], 'r--',color=tableau20[3],linewidth=3.0)
ax3.plot(data5[1], data5[2], 'r--',color=tableau20[5],linewidth=3.0)
ax3.plot(data6[1], data6[2], 'r--',color=tableau20[9],linewidth=3.0,label='KPR')
#ax3.scatter(out[1][0:n_out], out[2][0:n_out],color = 'r', marker = 'h',s = 40, label=' Hopfield')
#ax3.set_xlim([10**0.5,10**5.5])
ax3.set_ylim([10**(-4),10**2])
ax3.set_xlim([10**(1),10**6])
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xlabel('Dissipation$(k_B T/s)$', fontsize=25)
ax3.set_ylabel('Error rate', fontsize=25)


f.tight_layout()
plt.subplots_adjust(hspace = .1,wspace=0.2)
#f.subplots_adjust(left=None, bottom=None, right=None, top=None,
  #                  wspace=0.4, hspace=None)
f.set_size_inches(35, 10)
f.savefig('figure_T_adaption.png', bbox_inches='tight')






