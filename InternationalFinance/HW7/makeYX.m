% Making Y,X variables
function [Y0,YLm] = makeYX(Y,p)
  [T,k] = size(Y);
  Y0 = Y(p+1:T,:); % response variable

  % predictor variables which are the history of Y itself
  YL = zeros(T-p,p*k);
  for i = 1:p
    YL(:,k*(i-1)+1:k*i) = Y(p+1-i:T-i,:);
  end
  ki = p*k;
  kki = k*ki;
  YLm = zeros(k,kki,T-p);
  for t = 1:(T-p)
    YLm(:,:,t) = kron(eye(k),YL(t,:)); % p by k
  end
end