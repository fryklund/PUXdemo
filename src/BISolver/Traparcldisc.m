function [z,W,L] = Traparcldisc(curve,N,NL)
%[z,W] = Traparcldisc(f,N)
%[z,W] = Traparcldisc(f,N,NL)
%
%Discretizes the curve given by f(t) t = [0,2*pi] using the trapezoidal rule
%equidistant in arc-length .
%N is the number of quadrature points requested and if specified,
%NL is the number of 16-point Gauss-Legendre panels used to compute the 
%length of the curve. Default is NL = 1000.
%Non-adaptive. Make sure that the curve is well resolved by N points.

if nargin < 5
    NL = 1000;
end

[T,W] = GLinterval(0,2*pi,NL);
zp = curve.dtau(T);
L = W'*abs(zp);

L2 = linspace(1/N,1-1/N,N-1)'*L;
%Initial guess
t = linspace(0,2*pi,N+1)';
t = t(2:end-1);

dt = 1;
iter = 0;
while norm(dt)/norm(t) > 1e-13 && iter < 30
    %Set up 16-point GL quadrature on each segment.
    [T,W] = GLinterval(0,2*pi,N-1,[0;t]);
    zp = curve.dtau([t;T]);
    %Compute the cumulative sum of all the segment lengths
    F = cumsum(sum(reshape(W.*abs(zp(N:end)),16,N-1))');
    dt = (F-L2)./abs(zp(1:N-1));
    %Sort the parameters just in case. 
    t = sort(t - dt);
    
    iter = iter + 1;
end
if iter == 30
    warning('Traparcldisc : Newton did not converge, bad discretization');
end
t = [0;t];
z = curve.tau(t);
W = 2*pi/N*ones(N,1);
