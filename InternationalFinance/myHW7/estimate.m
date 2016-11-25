clc;
p    = 4;
mlag = 24;

Data    = xlsread('watson.xls',1,'A1:C167');

Y       = centre(Data);
m       = size(Y,2);
% set up data
X = zeros(T,m*p);
for t = 1:T
  X(t,:) = 


burnin  = 2000;
nSample = 2000;
thinin  = 4;

b     = 10;
D     = invpd(0.2*eye(m))/b;
D_inv = invpd(D);

mu_alpha        = zeros(m*m*p,1);
Sigma_alpha     = 1*eye(m*m*p);
% Compute constants before looping
Sigma_alpha_inv = inv(Sigma_alpha);
Sigma_inv_mu    = Sigma_alpha_inv*mu_alpha;


% initial values
Phi_old   = reshape(mu_alpha,m,p*m);
Omega_inv = b*D;


for i = 1:burnin
  [Phi,Fm,alpha]    = gen_Phi(X,y,Omega_inv,mu_alpha,Sigma_alpha_inv,Sigma_inv_mu,p,Phi_old);
  [Omega,Omega_inv] = gen_Omega(Y,X,beta,b,D_inv,Phi);
