function [u, uHat, fHat] = freeSpacePoissonFourierCoeff(N,L,f)
%% Solves u_xx + u_yy = f on [-L/2;L/2]^2

%% Initialisation
d = 2; % Number of spatial dimensions
sf = ceil(1+sqrt(d));%sampling factor %5/2
R = 1.5*L; %Size of truncation box: RxR
M = sf*N;
mid = M/2 + 1; % For clarity

disp(['Sampled points for free space Poisson: ' num2str(M) ])

%Define the truncated spectralrepresentation of Green's function
G0 = @(s2)  (1-besselj(0,R*sqrt(s2)))./s2 - R*log(R)*besselj(1,R*sqrt(s2))./sqrt(s2);
G0Limit = (1-2*log(R))*R^2/4; % From lim s--> 0 G0(s)

kMax = 2*pi*N/L/2; % Greatest wave represented on a N-points grid.
stepSize = 2*kMax/M;

%Create a new grid: MxM from kMax to kMax-stepSize. kMax is excluded due to
%periodicity

kv = -kMax + (0:M-1)/M*2*kMax;


[kX,kY] = meshgrid(kv,kv);


s2 = kX.^2 + kY.^2;

%% Create mollified Greens function

Ghatk = G0(s2);
Ghatk(mid,mid) = G0Limit;

%Shift to go from natural to swaped ordering. Then ifft, then go back to
%natural ordering.

G = real(ifftshift(ifft2(fftshift(Ghatk))));

%Pick out the 2*N points in the middle.
Gtruncated = G(mid-N:mid-1+N,mid-N:mid-1+N);

%Fourier transform the mollified Green's function
GhatkMoll = fft2(Gtruncated);

%% Solve Poisson with Mollified Greens

%Reshape f to a  NxN-matrix. Observe that it has natural ordering.
f = reshape(f,N,N);

%FFT and upsample by 2, i.e. pad with zeros. 
%Observe: no need to shift, since it will multiply with GhatkMoll, which has
%the same ordering.

fhatLong = fft2(f,2*N,2*N);

uHatLong = fhatLong.*GhatkMoll;


uLong = real(ifft2(uHatLong));
disp(['Imag(uLong) = ' num2str(norm(imag(uLong(N+1:2*N,N+1:2*N)))) ])
u = real(uLong(N+1:2*N,N+1:2*N));

%%%%%%%%%%%
% Remark: If you change which uHat for output, change shifting for boundary
% points for NUFFT.
% uHat = fft2((u)) should not be used. Not necessarily periodic. 
%%%%%%%%%%%

% uHatLong are the correct coefficients.
uHat = uHatLong;
fHat = fftshift(fft2(ifftshift(f)));


%save('/afs/kth.se/home/f/f/ffry/Private/research/poissonFuncExten/matlab/data/spectraPoisson.mat','fHat','uHat')


