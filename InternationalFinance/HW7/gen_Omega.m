function [Omega,Omega_inv] = gen_Omega(Y,X,beta,nu,R0)
  [T,k] = size(Y);
  ehat2 = zeros(k,k);

  for t = 1:T
    Xt = X(:,:,t);
    ehat = Y(t,:)' - Xt*beta;
    ehat2 = ehat2+ehat*ehat';
  end
  Omega1_inv = ehat2+invpd(R0);
  Omega1 = invpd(Omega1_inv);
  Omega_inv = randwishart(Omega1,(T+nu));
  Omega = invpd(Omega_inv);
end