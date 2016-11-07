function [ImpulseRespm,MHm] = recursive_VAR(n0,n1,Spec)
  b_   = Spec.b_;
  var_ = Spec.var_;
  p    = Spec.p;
  nu   = Spec.nu;
  R0   = Spec.R0;
  Y    = Spec.Y;
  mlag = Spec.mlag;

  k    = size(Y,2);


  Phi          = reshape(b_,p*k,k);
  Omega_inv    = nu*R0;
  ImpulseRespm = zeros(n1,mlag+1,k^2);

  pkk    = p*k*k;
  betam  = zeros(n1,pkk);
  Omegam = zeros(n1,k^2);

  [Y0,YLm] = makeYX(Y,p);
  n        = n0+n1;

  for iter = 1:n
    [~, resid] = minresid(iter,100);
    if resid == 0
      clc
      disp(['현재 반복시행은 ',num2str(iter)]);
    end

    [Phi,Fm,beta]     = gen_Phi(Y0,YLm,Phi,p,b_,var_,Omega_inv);
    [Omega,Omega_inv] = gen_Omega(Y0,YLm,beta,nu,R0);

    if iter > n0
      ImpulseRespm = gen_ImRes(Omega,Fm,mlag,n0,ImpulseRespm,iter);
      betam(iter-n0,:) = beta';
      Omegam(iter-n0,:) = vec(Omega)';
    end
  end
  MHm = [betam Omegam];

  plot_IRF(ImpulseRespm)
end
