function [Phi,Fm,beta] = gen_Phi(Y0,YLm,Phi0,p,b_,var_,Omega_inv)
  X = YLm;
  XX = 0;
  XY = 0;
  [T0,k] = size(Y0);
  for t = 1:T0
    Xt = X(:,:,t);
    XX = XX + Xt'*Omega_inv*Xt;
    XY = XY + Xt'*Omega_inv*Y0(t,:)';
  end
  precb_ = invpd(var_);
  B1_inv = precb_ + XX;
  B1_inv = 0.5*(B1_inv+B1_inv');
  B1 = invpd(B1_inv);
  B1 = 0.5*(B1+B1');
  A = XY + precb_*b_;
  BA = B1*A;

  chol_B1 = cholmod(B1)';
  beta = BA+chol_B1*randn(p*k*k,1);

  Phi = reshape(beta,p*k,k);
  Fm = [Phi'; eye((p-1)*k), zeros(k*(p-1),k)];

  eigF = eig(Fm);
  if maxc(abs(eigF)) >= 1
    Phi = Phi0
    Fm = [Phi'; eye((p-1)*k), zeros(k*(p-1),k)];
  end
end
