function plot_IRF(impulse_response)
  [~,mlag1,k2] = size(impulse_response); % k2 = k^2, mlag1 = mlag + 1
  k = sqrt(k2);
  ql = [0.025; 0.5; 0.975]; %For 5% credible interval
  xa = 0:(mlag1-1);

  a = reshape(1:k2,k,k);
  figure
  zeroline = zeros(mlag1,1);
  for i = 1:k2
    % impulse_response_ij = impulse_response(:,:,i); % n1 by (mlag+1)
    impulse_response_ij = quantile(impulse_response(:,:,i),ql)'; % (mlag+1) by 3

    [r,c] = find(a==i);
    subplot(k,k,i);
    plot(xa,impulse_response_ij(:,1),'k--',xa,impulse_response_ij(:,2),'b-',xa,impulse_response_ij(:,3),'k--',xa,zeroline,'k:','linewidth',2);
    xlim([0,mlag1]);
    title(['shock ', num2str(c), ' to vari ', num2str(r)])
  end
end
