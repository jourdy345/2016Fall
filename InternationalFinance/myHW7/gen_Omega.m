function [Omega,Omega_inv] = gen_Omega(Y,X,beta,b,D_inv,Phi)
  [T,m] = size(Y);
  
  D_T = invpd(D_inv+(Y'-Phi*X')*(Y-X*B'));
  Omega_inv = randwishart(D_T,(T+b));
  Omega = invpd(Omega_inv);
  
end