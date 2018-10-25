function [bc, f] = example3()

bc = @(p) exp(-1)*(cos(p(:,1))+cos(p(:,2))) + exp(-sqrt(2))*(cos(p(:,1)*2)+cos(p(:,2)*2)) + exp(-sqrt(4))*(cos(p(:,1)*4)+cos(p(:,2)*4)) + exp(-sqrt(8))*(cos(p(:,1)*8)+cos(p(:,2)*8)) + exp(-sqrt(16))*(cos(p(:,1)*16)+cos(p(:,2)*16)) + exp(-sqrt(32))*(cos(p(:,1)*32)+cos(p(:,2)*32));
% f = \Delta bc.
f = @(p) -exp(-1)*(cos(p(:,1))+cos(p(:,2))) - exp(-sqrt(2))*(4*cos(p(:,1)*2)+4*cos(p(:,2)*2)) - exp(-sqrt(4))*(16*cos(p(:,1)*4)+16*cos(p(:,2)*4)) - exp(-sqrt(8))*(64*cos(p(:,1)*8)+64*cos(p(:,2)*8)) - exp(-sqrt(16))*(256*cos(p(:,1)*16)+256*cos(p(:,2)*16)) - exp(-sqrt(32))*(1024*cos(p(:,1)*32)+1024*cos(p(:,2)*32));


end
