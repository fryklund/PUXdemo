function [bc,f] = example1()
bc = @(z) sin(2*pi*z(:,1)).*sin(2*pi*z(:,2))/(pi*pi*8);

% f = \Delta bc.
f = @(z) -sin(2*pi*z(:,1)).*sin(2*pi*z(:,2));
end