function [y] = centre(x)
  y = x-ones(size(x,1),1)*mean(x);
end