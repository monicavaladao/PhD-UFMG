function [f] = rosen(x)
d = length(x);
f=sum(100.*(x(:,1:d-1).^2-x(:,2:d)).^2+(x(:,1:d-1)-1).^2,2);

end