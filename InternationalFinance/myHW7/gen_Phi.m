function [Phi,Fm,alpha] = gen_Phi(X,y,Omega_inv,mu_alpha,Sigma_alpha_inv,Sigma_inv_mu,p,Phi_old)
  Xcol = size(X,2);
  [~,m] = size(Omega_inv);
  Sigma_T = invpd(kron(Omega_inv,X'*X)+Sigma_alpha_inv);
  mu_T    = Sigma_T*(Sigma_inv_mu+kron(Omega_inv,X')*y);
  alpha = mu_T+chol(Sigma_T,'lower')*randn(m*Xcol,1);
  Phi = reshape(alpha,m,Xcol);
  Fm = [Phi; eye((p-1)*m),zeros((p-1)*m,m)];
  eigF = eig(Fm);
  if maxc(abs(eigF)) > = 1
    Phi = Phi_old;
    Fm = [Phi; eye((p-1)*m),zeros((p-1)*m,m)];
  end
end

%   X = YLm;
%   XX = 0;
%   XY = 0;
%   [T0,k] = size(Y0);
%   for t = 1:T0
%     Xt = X(:,:,t);
%     XX = XX + Xt'*Omega_inv*Xt;
%     XY = XY + Xt'*Omega_inv*Y0(t,:)';
%   end
%   precb_ = invpd(var_);
%   B1_inv = precb_ + XX;
%   B1_inv = 0.5*(B1_inv+B1_inv');
%   B1 = invpd(B1_inv);
%   B1 = 0.5*(B1+B1');
%   A = XY + precb_*b_;
%   BA = B1*A;

%   chol_B1 = cholmod(B1)';
%   beta = BA+chol_B1*randn(p*k*k,1);

%   Phi = reshape(beta,p*k,k);
%   Fm = [Phi'; eye((p-1)*k), zeros(k*(p-1),k)];

%   eigF = eig(Fm);
%   if maxc(abs(eigF)) >= 1
%     Phi = Phi0
%     Fm = [Phi'; eye((p-1)*k), zeros(k*(p-1),k)];
%   end
% end
