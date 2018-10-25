function [muLapl,ALapl] = dlpDensitySolver(N,M,idxBody,z,dz,d2,w,RHS,S)
% Calculate density mu from the boundary integral formulation for Laplace
% eq. with Dirichlet boundary condition (RHS)
% Singular integrals calcualted with limit points
%Notation follows from http://www.sciencedirect.com/science/article/pii/S0021999183710739?via%3Dihub

K = eye(N,N);
B = zeros(N,M);
C = zeros(M,N);
D = zeros(M,M);
for i=1:N
    K(i,:) = K(i,:) + (1/pi*w.*imag(dz./(z-z(i))))';
end
d = 1 + 1/pi*w.*imag(d2./(2*dz));
b = [2*RHS;zeros(M,1)];%imag(tauk.^2./(tauk-zp));
K(logical(eye(size(K)))) = 0;
K = K + diag(d,0);

for i = 1:M
B(:,i) = B(:,i) + 2 * log(abs(z - S(i)));
end

for i = 1:M
    idx = (1 + idxBody(i+1)):idxBody(i+2);
    temp_row = zeros(1,N);
    temp_row(idx) = 2 * w(idx).*abs(dz(idx));
   C(i,:) = C(i,:) + temp_row;
end

A = [K,B;C,D];

muLapl = gmres(A,b,[],1e-14,N);
ALapl = muLapl((N+1):end);
muLapl = muLapl(1:N);
disp(['Size of coefficient: ' num2str(abs(ALapl))])
disp(['Integral of mu over interior boundary curves: ' num2str(sum(muLapl((1+idxBody(2)):end).*2 .* w((1+idxBody(2)):end).*abs(dz((1+idxBody(2)):end))))]);
% mu_lapl = A\b;
end