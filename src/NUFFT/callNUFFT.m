function [u] = callNUFFT(uHat,M,eval_pnts,LBox,tol)

iflag = +1;
os = size(uHat,1)/M;%Oversampling
pointsShift = [real(eval_pnts)'*pi/LBox-pi;imag(eval_pnts)'*pi/LBox-pi]/os;
% .' is transpose, not hermitian transpose '

uNufft = (uHat/(os*M)^2).';
u = real(nufft2d2(size(pointsShift,2),pointsShift(1,:),pointsShift(2,:),iflag,tol,os*M,os*M,uNufft));

end