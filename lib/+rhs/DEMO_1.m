function [bc,f] = DEMO_1()
bc = @(z)  -exp(-(sqrt(z(:,1).^2 + z(:,2).^2)).^2);

% f = \Delta bc.
f = @(z) exp(-z(:,1).^2-z(:,2).^2) * 4 .*(-z(:,1).^2 - z(:,2).^2 + 1);



end