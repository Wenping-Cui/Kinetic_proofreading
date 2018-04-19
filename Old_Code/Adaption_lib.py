# Functions to calculate energy, speed and error rate
def model(N,m,para):
# This is a specific case for N=4, m=2
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
#    Master_f = np.array([[-kappa*R-(1.0/tau+phi),b-kappa*R,-kappa*R,-kappa*R,-kappa*R],[phi-gamma*R,-gamma*R-(1.0/tau+phi+b),b-gamma*R,-gamma*R,-gamma*R],[-gamma*R,phi-gamma*R,-gamma*R-(1.0/tau+phi+b),b-gamma*R,-gamma*R],[-gamma*R,-gamma*R,phi-gamma*R,-gamma*R-(1.0/tau+alpha*K+b),b-gamma*R], [-gamma*R,-gamma*R,-gamma*R,alpha*K-gamma*R,-gamma*R-(b+1.0/tau)]]);
#     B_f = np.array([-kappa*R*L,-gamma*R*L,-gamma*R*L,-gamma*R*L,-gamma*R*L]);
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
    h = 10**(-10)*phi; # Step size
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
       Energy_C = (phi*C[0]-b*C[1])*np.log(kappa*tauf)+J_C*np.log(phi/b/gamma)+(alpha*K*C[N-1]-b*C[N])*np.log(alpha*K/b/gamma);
# For D
    J_D = 0;
    for j in np.arange(N-1):
        J_D = J_D + phi*D[j]-b*D[j+1];
    Energy_D = 0;
    if gamma!=0:
       Energy_D = (phi*D[0]-b*D[1])*np.log(kappa*taus)+J_D*np.log(phi/b/gamma)+(alpha*K*D[N-1]-b*D[N])*np.log(alpha*K/b/gamma);
# Total energy consumption
    energy = Energy_C+Energy_D;
    TJ = 1.0/(alpha*K*C[N-1]-b*C[N])
    return  T, energy, eta, D[N]+C[N]
#    return  eta, J0_D-J1_D, J1_D-J2_D, J2_D-J3_D,J0_D,J1_D,J2_D,J3_D
#    return  D[0], D[1], D[2], D[4],C[0], C[1], C[2], C[4] 

def model_simple(N,m,para):
# This is a specific case for N=4, m=2
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
          Master_f = np.array([[-kappa*R-(1.0/tau+phi),b-kappa*R,-kappa*R,-kappa*R,-kappa*R],[phi-gamma*R/tauf,-gamma*R/tauf-(1.0/tau+phi+b),b-gamma*R/tauf,-gamma*R/tauf,-gamma*R/tauf],[-gamma**2*R/tauf,phi-gamma**2*R/tauf,-gamma**2*R/tauf-(1.0/tau+phi+b),b-gamma**2*R/tauf,-gamma**2*R/tauf],[-gamma**3*R/tauf,-gamma**3*R/tauf,phi-gamma**3*R/tauf,-gamma**3*R/tauf-(1.0/tau+alpha*K+b),b-gamma**3*R/tauf], [-gamma**4*R/tauf,-gamma**4*R/tauf,-gamma**4*R/tauf,alpha*K-gamma**4*R/tauf,-gamma**4*R/tauf-(b+1.0/tau)]]);
          B_f = np.array([-kappa*R*L,-gamma*R*L/tauf,-gamma**2*R*L/tauf,-gamma**3*R*L/tauf,-gamma**4*R*L/tauf]);
          C = np.linalg.solve(Master_f, B_f)
# Solve the master equations for D wrong product
          tau = taus;
          L = Ls;
          Master_s = np.array([[-kappa*R-(1.0/tau+phi),b-kappa*R,-kappa*R,-kappa*R,-kappa*R],[phi-gamma*R/tau,-gamma*R/tau-(1.0/tau+phi+b),b-gamma*R/tau,-gamma*R/tau,-gamma*R/tau],[-gamma**2*R/tau,phi-gamma**2*R/tau,-gamma**2*R/tau-(1.0/tau+phi+b),b-gamma**2*R/tau,-gamma**2*R/tau],[-gamma**3*R/tau,-gamma**3*R/tau,phi-gamma**3*R/tau,-gamma**3*R/tau-(1.0/tau+alpha*K+b),b-gamma**3*R/tau], [-gamma**4*R/tau,-gamma**4*R/tau,-gamma**4*R/tau,alpha*K-gamma**4*R/tau,-gamma**4*R/tau-(b+1.0/tau)]]);
          B_s = np.array([-kappa*R*L,-gamma*R*L/tau,-gamma**2*R*L/tau,-gamma**3*R*L/tau,-gamma**4*R*L/tau]);
          D = np.linalg.solve(Master_s, B_s)
# Check the code and it works at b=0
#    error_acc = C[N]/C[0] - epsilon/sigma*alpha*KT*phi*tauf**2*(phi*tauf/(phi*tauf+1.0))**2/(epsilon/sigma*(1.0+alpha*KT*tauf)+(phi*tauf/(phi*tauf+1.0))**2*C[0]+(phi*taus/(phi*taus+1.0))**2*D[0]);
# Calculate the error rate: wrong product/right product
    eta = D[4]/C[4];


# Calculate the first passage time
    h = 10**(-5)*phi; # Step size
    K = KT/(1.0 + sigma/epsilon*(C[m]+D[m]));
    tau = taus;
    L = Ls;
    A_s = np.array([[-kappa*R-(gamma/tau+gamma**2/tau+gamma**3/tau+gamma**4/tau)*R,1.0/tau,1.0/tau,1.0/tau,1.0/tau,1.0/tau,0],[kappa*R,-1.0/tau-phi,b,0,0,0,0],[gamma*R/tau,phi,-1.0/tau-phi-b,b,0,0,0],[gamma**2*R/tau,0,phi,-1.0/tau-phi-b,b,0,0],[gamma**3*R/tau,0,0,phi+alpha*K*sigma*C[3]/(epsilon+sigma*(C[2]+D[2])),-1.0/tau-alpha*K-b,b,0],[gamma**4*R/tau,0,0,-alpha*K*sigma*C[3]/(epsilon+sigma*(C[2]+D[2])),alpha*K,-b-1.0/tau-W,0],[0,0,0,0,0,W,0]]);
    B_s = h*np.eye(A_s.shape[0]);   
    C_s = np.zeros(A_s.shape[0]);
    C_s[0] = 1.0;

    A_s1 = np.array([[-kappa*R-(gamma+gamma**2+gamma**3+gamma**4)*R/tau,1.0/tau,1.0/tau,1.0/tau,1.0/tau,1.0/tau],[kappa*R,-1.0/tau-phi,b,0,0,0],[gamma*R/tau,phi,-1.0/tau-phi-b,b,0,0],[gamma**2*R/tau,0,phi,-1.0/tau-phi-b,b,0],[gamma**3*R/tau,0,0,phi+alpha*K*sigma*C[3]/(epsilon+sigma*(C[2]+D[2])),-1.0/tau-alpha*K-b,b],[gamma**4*R/tau,0,0,-alpha*K*sigma*C[3]/(epsilon+sigma*(C[2]+D[2])),alpha*K,-b-1.0/tau-W]]);
    C_s1 = np.zeros(A_s.shape[0]-1);
    C_s1[0] = 1.0;
# The first passge time in the unit of 1/W, T should be independent of W
    T_s = -W**2*(np.linalg.solve(B_s-A_s, C_s)[A_s.shape[0]-2]-np.linalg.solve(-A_s1, C_s1)[A_s.shape[0]-2])/h;
    speed = T_s;


# Calculate the first passage time
    h = 2*10**(-8)*phi; # Step size
    K = KT/(1.0 + sigma/epsilon*(C[m]+D[m]));
    tau = tauf;
    L = Lf;
    A_f = np.array([[-kappa*R-(gamma/tau+gamma**2/tau+gamma**3/tau+gamma**4/tau)*R,1.0/tau,1.0/tau,1.0/tau,1.0/tau,1.0/tau,0],[kappa*R,-1.0/tau-phi,b,0,0,0,0],[gamma*R/tau,phi,-1.0/tau-phi-b,b,0,0,0],[gamma**2*R/tau,0,phi,-1.0/tau-phi-b,b,0,0],[gamma**3*R/tau,0,0,phi+alpha*K*sigma*C[3]/(epsilon+sigma*(C[2]+D[2])),-1.0/tau-alpha*K-b,b,0],[gamma**4*R/tau,0,0,-alpha*K*sigma*C[3]/(epsilon+sigma*(C[2]+D[2])),alpha*K,-b-1.0/tau-W,0],[0,0,0,0,0,W,0]]);
    B_f = h*np.eye(A_f.shape[0]);   
    C_f = np.zeros(A_f.shape[0]);
    C_f[0] = 1.0;

    A_f1 = np.array([[-kappa*R-(gamma+gamma**2+gamma**3+gamma**4)*R/tau,1.0/tau,1.0/tau,1.0/tau,1.0/tau,1.0/tau],[kappa*R,-1.0/tau-phi,b,0,0,0],[gamma*R/tau,phi,-1.0/tau-phi-b,b,0,0],[gamma**2*R/tau,0,phi,-1.0/tau-phi-b,b,0],[gamma**3*R/tau,0,0,phi+alpha*K*sigma*C[3]/(epsilon+sigma*(C[2]+D[2])),-1.0/tau-alpha*K-b,b],[gamma**4*R/tau,0,0,-alpha*K*sigma*C[3]/(epsilon+sigma*(C[2]+D[2])),alpha*K,-b-1.0/tau-W]]);
    C_f1 = np.zeros(A_f.shape[0]-1);
    C_f1[0] = 1.0;
# The first passge time in the unit of 1/W, T should be independent of W
    T_f = -W**2*(np.linalg.solve(B_f-A_f, C_f)[A_f.shape[0]-2]-np.linalg.solve(-A_f1, C_f1)[A_f.shape[0]-2])/h;
    speed = T_f;
# Calculate the energy consumption
# For C
    K = KT/(1.0 + sigma/epsilon*(C[m]+D[m]));
    L_C = L - C[0] -C[1]-C[2]-C[3]-C[4]
    J0_C = phi*C[0]-b*C[1];
    J1_C = phi*C[1]-b*C[2];
    J2_C = phi*C[2]-b*C[3];
    J3_C = alpha*K*C[3]-b*C[4];
    Energy_C =0;
    if gamma !=0:
       Energy_C = J0_C*np.log(kappa*phi/gamma/b)+J1_C*np.log(phi/(b*gamma))+J2_C*np.log(phi/(b*gamma))+J3_C*np.log(alpha*K/(b*gamma));
# For D
    L_D = L - D[0] -D[1]-D[2]-D[3]-D[4]
    J0_D = phi*D[0]-b*D[1];
    J1_D = phi*D[1]-b*D[2];
    J2_D = phi*D[2]-b*D[3];
    J3_D = alpha*K*D[3]-b*D[4];
    Energy_D = 0;
    if gamma !=0:
       Energy_D = J0_D*np.log(kappa*phi/gamma/b)+J1_D*np.log(phi/(b*gamma))+J2_D*np.log(phi/(b*gamma))+J3_D*np.log(phi/(b*gamma));
#    Energy_D = J0_D*np.log(kappa/gamma)+(J0_D+J1_D+J2_D)*np.log(phi/b)+J3_D*np.log(alpha*K/b);
# Total energy consumption
    energy = Energy_C+Energy_D;
    return  T_f, energy, eta,T_s, C, D
#    return  eta, J0_D-J1_D, J1_D-J2_D, J2_D-J3_D,J0_D,J1_D,J2_D,J3_D
#    return  D[0], D[1], D[2], D[4],C[0], C[1], C[2], C[4] 

