# Functions to calculate energy, speed and error rate
def model(N,m,para):
    import numpy as np
# para = np.array([Ls,Lf,taus,tauf,R,kappa,KT,sigma,epsilon,alpha,b,phi,gamma,W])
    C = np.zeros(N+1);  # Agonist complex phosphorylated n times
    D = np.zeros(N+1);  # Non-Agonist complex phosphorylated n times

    Ls = para[0];    # Concentration of self ligands 10^4 - 10^5
    Lf = para[1];       # Concentration of foreign ligands ~ 10 per cell
    taus = para[2];      # Self complex bonding time
    tauf = para[3];     # Foreign complex bonding time

    R = para[4];     # Concentration of receptor when Ls = 0
    kappa = para[5]; # Receptor-ligand bonding rate

    KT = para[6];    # Kinase total 
    sigma = para[7];   # Kinase phosphorylation rate
    epsilon = para[8]; # Kinase dephosphorylation rate
    K = KT/(1.0 + sigma/epsilon*(C[m]+D[m]));  # Kinase

    alpha = para[9]; # alpha*K  the last forward rate 
    b = para[10];  # Back rate
    phi = para[11]; # Forward rate
    gamma = para[12]; 
    W = para[13]; # Reaction rate to the final product
    # Solve the master equations with iterative method
    i=1.0;
    while np.abs(i)<10:
          i = i+1;
# K works as a coupling between foreign ligand C and self ligand D
          K = KT/(1.0 + sigma/epsilon*(C[m]+D[m]));
# Solve the master equations for C correct product
          tau = tauf;
          L = Lf;
          Master_f = np.zeros((N+1,N+1));
          for j in np.arange(N+1):
              Master_f[j,:] = -gamma**j*R/tau;
          Master_f[0,:] = -kappa*R;
          for j in np.arange(N):
              Master_f[j, j+1] = Master_f[j, j+1] + b;
              Master_f[j+1, j] = Master_f[j+1, j] + phi;
              Master_f[j,j] = Master_f[j,j] - (1.0/tau + phi +b);
          Master_f[0,0] = -kappa*R - (1.0/tau + phi);
          Master_f[N,N] = -gamma**N*R/tau - (1.0/tau + b);
          Master_f[N-1,N-1] = -gamma**(N-1)*R/tau -(1.0/tau + b+ alpha*K);
          Master_f[N,N-1] = -gamma**N*R/tau + alpha*K;
          B_f = np.zeros(N+1);
          for j in np.arange(N):
              B_f[j+1] = -gamma**(j+1)*R*L/tau;
          B_f[0] = -kappa*R*L;
          C = np.linalg.solve(Master_f, B_f)
# Solve the master equations for D wrong product
          tau = taus;
          L = Ls;
          Master_s = np.zeros((N+1,N+1));
          for j in np.arange(N+1):
              Master_s[j,:] = -gamma**j*R/tau;
          Master_s[0,:] = -kappa*R;
          for j in np.arange(N):
               Master_s[j, j+1] = Master_s[j, j+1] + b;
               Master_s[j+1, j] = Master_s[j+1, j] + phi;
               Master_s[j,j] = Master_s[j,j] - (1.0/tau + phi +b);
          Master_s[0,0] = -kappa*R - (1.0/tau + phi);
          Master_s[N,N] = -gamma**N*R/tau - (1.0/tau + b);
          Master_s[N-1,N-1] = -gamma**(N-1)*R/tau -(1.0/tau + b+ alpha*K);
          Master_s[N,N-1] = -gamma**N*R/tau + alpha*K;
          B_s = np.zeros(N+1);
          for j in np.arange(N):
              B_s[j+1] = -gamma**(j+1)*R*L/tau;
          B_s[0] = -kappa*R*L;
          D = np.linalg.solve(Master_s, B_s)
# Check the code and it works at b=0
#    error_acc = C[N]/C[0] - epsilon/sigma*alpha*KT*phi*tauf**2*(phi*tauf/(phi*tauf+1.0))**(N-2)/(epsilon/sigma*(1.0+alpha*KT*tauf)+(phi*tauf/(phi*tauf+1.0))**m*C[0]+(phi*taus/(phi*taus+1.0))**m*D[0]);


# Calculate the error rate: wrong product/right product
    eta = D[N]/C[N];


# Calculate the first passage time
    h = 10**(-8)*phi; # Step size
    tau = tauf;
    K = KT/(1.0 + sigma/epsilon*(C[m]+D[m]));
    A_s = np.zeros((N+1+2,N+1+2))
    A_s[0,:] = 1.0/tauf;
    A_s[0,N+2] = 0;
    A_s[1,0] = kappa*R
    A_s[0,0] = -kappa*R
    for j in np.arange(N):
        A_s[0,0] = A_s[0,0] -gamma**(j+1)*R/tau 
        A_s[j+2,0] = gamma**(j+1)*R/tau
    A_s[N+2,0] =0;
    for j in np.arange(1,N+1):
              A_s[j, j+1] = b;
              A_s[j+1, j] = phi;
              A_s[j,j] = - (1.0/tauf + phi +b);
    A_s[1,1] = -(1.0/tauf + phi);
    A_s[N+1,N+1] = - (1.0/tauf +b +W);
    A_s[N,N] = -(1.0/tauf+alpha*K+b);
    A_s[N+1,N] = alpha*K;        
    A_s[N+2,N+1] = W;
    A_s[N+1,m+1] = A_s[N+1,m+1]- alpha*K*sigma*C[N-1]/(epsilon+(sigma*C[m]+sigma*D[m]));
    A_s[N,m+1] = A_s[N,m+1] + alpha*K*sigma*C[N-1]/(epsilon+sigma*(C[m]+D[m])) ;
    B_s = h*np.eye(A_s.shape[0]);   
    C_s = np.zeros(A_s.shape[0]);
    C_s[0] = 1.0;

    A_s1 = np.zeros((N+1+1,N+1+1))
    A_s1[0,:] = 1.0/tauf;
    A_s1[1,0] = kappa*R
    A_s1[0,0] = -kappa*R
    for j in np.arange(N):
        A_s1[0,0] = A_s1[0,0] -gamma**(j+1)*R/tau 
        A_s1[j+2,0] = gamma**(j+1)*R/tau 
    for j in np.arange(1,N+1):
              A_s1[j, j+1] = b;
              A_s1[j+1, j] = phi;
              A_s1[j,j] = - (1.0/tauf + phi +b);
    A_s1[1,1] = -(1.0/tauf + phi);
    A_s1[N+1,N+1] = - (1.0/tauf +b +W);
    A_s1[N,N] = -(1.0/tauf+alpha*K+b);
    A_s1[N+1,N] = alpha*K;        
    A_s1[N+1,m+1] = A_s1[N+1,m+1]- alpha*K*sigma*C[N-1]/(epsilon+(sigma*C[m]+sigma*D[m]));
    A_s1[N,m+1] = A_s1[N,m+1] + alpha*K*sigma*C[N-1]/(epsilon+sigma*(C[m]+D[m]));
    C_s1 = np.zeros(A_s.shape[0]-1);
    C_s1[0] = 1.0;
# The first passge time in the unit of 1/W, T should be independent of W
    T = -W*(np.linalg.solve(B_s-A_s, C_s)[A_s.shape[0]-2]-np.linalg.solve(-A_s1, C_s1)[A_s.shape[0]-2])/h;


# Calculate the energy consumption
# For C
    K = KT/(1.0 + sigma/epsilon*(C[m]+D[m]));
    J_C = 0;
    for j in np.arange(N-1):
        J_C = J_C + phi*C[j]-b*C[j+1];
    Energy_C = 0;
    if gamma!=0:
       Energy_C = (phi*C[0]-b*C[1])*np.log(kappa)+J_C*np.log(phi/b/gamma)+(alpha*K*C[N-1]-b*C[N])*np.log(alpha*K/b/gamma);
# For D
    J_D = 0;
    for j in np.arange(N-1):
        J_D = J_D + phi*D[j]-b*D[j+1];
    Energy_D = 0;
    if gamma!=0:
       Energy_D = (phi*D[0]-b*D[1])*np.log(kappa)+J_D*np.log(phi/b/gamma)+(alpha*K*D[N-1]-b*D[N])*np.log(alpha*K/b/gamma);
# Total energy consumption
    energy = Energy_C+Energy_D;
    return  T, energy, eta,T, C,D

