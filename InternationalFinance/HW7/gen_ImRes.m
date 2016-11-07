function ImpulseRespm = gen_ImRes(Omega,F,mlag,n0,ImpulseRespm,iter)
  Binv = chol(Omega)';
  k = size(Omega,1);
  FF = eye(size(F,1));
  for j = 1:(mlag+1)
    psi_j = FF(1:k,1:k);
    theta = psi_j*Binv;
    theta = vec(theta);
    for i = 1:k^2
      ImpulseRespm(iter-n0,j,i) = theta(i);
    end
    FF = FF*F;
  end
end
