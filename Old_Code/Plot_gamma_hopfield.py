# Code to solve parallel adaptive sorting 
# Specific case: N=4 and m=2
import numpy as np
import matplotlib.pyplot as plt
from Adaption_lib import *
from matplotlib import colors
import matplotlib.ticker as mticker
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

Ls = 10**4;    # Concentration of self ligands 10^4 - 10^5
Lf = 10**4;       # Concentration of foreign ligands ~ 10 per cell
taus = 1.0;      # Self complex bonding time
tauf = 10.0;     # Foreign complex bonding time

R = 10.0**4;     # Concentration of receptor when Ls = 0
kappa = 3.0*10**(2); # Receptor-ligand bonding rate
kappaR = 3.0;    # Considering Ls, suppose the receptors are largely in excess

KT = 1.0*10**(4);    # Kinase total
KT = 1.0*10**(4);    # Kinase total
sigma = 0;   # Kinase phosphorylation rate
epsilon = 0.01; # Kinase dephosphorylation rate
K = KT/(1.0 + sigma/epsilon*(C[m]+D[m]));  # Kinase


alpha = 3.0*10**(-4); # alpha*K  the last forward rate 
b = 0.1;  # Back rate
phi = 1.0; # Forward rate
gamma = 10**(-3); 
W = 0.01; # Reaction rate to the final product
# Summary the parameters

para = np.array([Ls,Lf,taus,tauf,R,kappa,KT,sigma,epsilon,alpha,b,phi,gamma,W])

# Output: speed, energy and eta
def data_ana(tau,gamma):
    Num = 150; # Number of scatter points
    data = np.zeros((Num,4+6));
    for num in range(0, Num):
#        kappa = 10**np.random.uniform(-4,4);    
#        sigma = 10**np.random.uniform(-2,0.5);
#        epsilon = 10**np.random.uniform(-2,0.5);
#        alpha = 10**np.random.uniform(-4,4);    
        phi = 10**(-4.6+(5.0+4.5)/Num*num);
        b = 0.01*phi;
#        kappa = phi/R;
#        sigma = 0;
        alpha = phi/KT
        para = np.array([Ls,Lf,tau,tauf,R,kappa,KT,sigma,epsilon,alpha,b,phi,gamma,W])
        output = model(N,m,para);
        data[num] = np.concatenate((output[0:4],[kappa,sigma,epsilon,alpha,b,phi]), axis=0)
    return data[:,0:4+6]

# Sorting by eta(error rate)

# Define bin size(of the array size) for energy consumption

data1 = data_ana(1.0,10**(-1));
data2 = data_ana(1.0,10**(-1.5));
data3 = data_ana(1.0,10**(-3));
#data4 = data_ana(0.001);
n_out = 1
#data = np.concatenate((data, data4), axis=0)
data1 = data1[data1[:,2].argsort()]
data2 = data2[data2[:,2].argsort()]
data3 = data3[data3[:,2].argsort()]
out1 = data1[0:n_out]
out1 = np.transpose(out1);
out2 = data2[0:n_out]
out2 = np.transpose(out2);
out3 = data3[0:n_out]
out3 = np.transpose(out3);
out = np.concatenate((data1[0:n_out], data2[0:n_out]), axis=0)
out = np.concatenate((out, data3[0:n_out]), axis=0)

out = np.transpose(out);
data1 = data_ana(1.0,10**(-1));
data2 = data_ana(1.0,10**(-1.5));
data3 = data_ana(1.0,10**(-3));
data4 = data_ana(1.0,0);
data1 = np.transpose(data1);
data2 = np.transpose(data2);
data3 = np.transpose(data3); 
data4 = np.transpose(data4); 
#data4 = np.transpose(data4);
bin_error = (taus/tauf)**(N) ;

print bin_error/out[2][0:n_out]
print bin_error



f, ((ax1, ax2, ax3)) = plt.subplots(1, 3, sharex=False, sharey=False)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 25})
ax1.plot(data1[0], data1[2],color=tableau20[2],linewidth=3.0, label=' $\\gamma = 10^{-1}$')
ax1.plot(data2[0], data2[2],color=tableau20[4],linewidth=3.0,label=' $\\gamma = 10^{-1.5}$')
ax1.plot(data3[0], data3[2],color=tableau20[8],linewidth=3.0, label=' $\\gamma = 10^{-3}$')
ax1.plot(data4[0], data4[2],color=tableau20[6], linestyle='-',dashes = [5,2,5,2],linewidth=3.0, label=' $\\gamma = 0$')
ax1.scatter(out1[0], out1[2],color=tableau20[2], marker = 's',s = 100)
ax1.scatter(out2[0], out2[2],color=tableau20[4], marker = 's',s = 100)
ax1.scatter(out3[0], out3[2],color=tableau20[8], marker = 's',s = 100)

#ax1.scatter(data4[0], data4[2],color = 'm', marker = '^', label=' $b = 0.001\\phi$')
ax1.set_xscale('log')
#ax1.scatter(speed_boundary, eta_boundary, s=20, c='r')
ax1.set_yscale('log')
#ax1.set_xlabel('Time', fontsize=25)
legend = ax1.legend(loc='upper left',frameon=False, fontsize=20)

ax1.set_xlabel('MFPT$(sec)$', fontsize=25)
ax1.set_ylabel('Error rate', fontsize=25)





ax2.plot(data1[0], data1[1],color=tableau20[2],linewidth=3.0, label=' $\\gamma = 10^{-1}$')
ax2.plot(data2[0], data2[1],color=tableau20[4],linewidth=3.0,label=' $\\gamma = 10^{-1.5}$')
ax2.plot(data3[0], data3[1],color=tableau20[8],linewidth=3.0, label=' $\\gamma = 10^{-3}$')
ax2.scatter(out1[0], out1[1],color=tableau20[2], marker = 's',s = 100)
ax2.scatter(out2[0], out2[1],color=tableau20[4], marker = 's',s = 100)
ax2.scatter(out3[0], out3[1],color=tableau20[8], marker = 's',s = 100)
ax2.set_ylim([10**1,10**5])
ax2.set_xscale('log')
ax2.set_yscale('log')
legend = ax2.legend(loc='upper right',frameon=False, fontsize=20)
ax2.set_xlabel('MFPT$(sec)$', fontsize=25)
ax2.set_ylabel('Dissipation$(k_B T/s)$', fontsize=25)

ax3.plot(data1[1], data1[2],color=tableau20[2],linewidth=3.0, label=' $\\gamma = 10^{-1}$')
ax3.plot(data2[1], data2[2],color=tableau20[4],linewidth=3.0,label=' $\\gamma = 10^{-1.5}$')
ax3.plot(data3[1], data3[2],color=tableau20[8],linewidth=3.0, label=' $\\gamma = 10^{-3}$')
ax3.scatter(out1[1], out1[2],color=tableau20[2], marker = 's',s = 100)
ax3.scatter(out2[1], out2[2],color=tableau20[4], marker = 's',s = 100)
ax3.scatter(out3[1], out3[2],color=tableau20[8], marker = 's',s = 100)
ax3.set_yscale('log')
ax3.set_xscale('log')

legend = ax3.legend(loc='upper right',frameon=False, fontsize=20)
ax3.set_xlabel('Dissipation$(k_B T)$', fontsize=25)
ax3.set_ylabel('Error rate', fontsize=25)



f.tight_layout()

plt.subplots_adjust(hspace = .1,wspace=0.4)
#f.subplots_adjust(left=None, bottom=None, right=None, top=None,
  #                  wspace=0.4, hspace=None)
f.set_size_inches(35, 10)
f.savefig('figure_gamma_hopfield.png', bbox_inches='tight')






