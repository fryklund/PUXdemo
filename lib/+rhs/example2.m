function [bc,f] = example2()
c = 10;
bc = @(z) sin(c*sum(z,2))+z(:,1).^2/9-z(:,2)+8+exp(-5*c^2*z(:,1).^2);

% f = \Delta bc.
f = @(z) -2*c^2*sin(c*sum(z,2))+2/9+...
    10*c^2*(10*c^2*z(:,1).^2-1).*exp(-5*c^2*z(:,1).^2);
end

