function [errorL2,errorMax] = example2(varargin)
%% The Poisson equation - Mutiply connected
% Solves the Poisson equation on a mutiply connected domain. See Example 2
% in https://arxiv.org/abs/1712.08461.
%
%    u - \alpha^2 \Delta u = f in Omega
%     u = g on \Gamma,
%
%

%
% Requires the following to run:
% RBF-QR by Courant Mathematics and Computing Laboratory. See Readme for links
% NUFFT by Elisabeth Larsson. See Readme for links.
%
%%

switch nargin
    case 0
    M = 800;
    params.regularity = 0;
    case 1
        params.M = varargin{1};
        params.regularity = 0;
    case 2
    params.M = varargin{1};
    params.regularity = varargin{2};
    otherwise
        warning('Error: either no input or resolution, RBF and nbr RBF-centers')
        return;
end

% Store parameters for PUX in PUXstruct. See setupPUX for explanation of 
% input.
PUXstruct.params = struct('ep',2,'coarsen',2,'R_ratio',2.5,'R_p',0.0675,'regularity',params.regularity);


%% Load stored data
% Load interpolatory quadrature weights for interpolating from 16 points
% panels to 32, and vice versa.
load 'IP1632.mat'

% Load Gaus-Legendre weigths for  16 and 32 points.
load 'glW.mat' 




%% Set up computational domain.
% Set size of box B = [-LBox,LBox]^2.
LBox = 0.4;

% Setup geometry. The function diffCurve gives symbolic expressions for 
%the first and second derivative of the curve, the normal, the tangent and
% the curvature. 

curve(1) = diffCurve(curves.example2());
curve(2) = diffCurve(curves.example2Exterior());

% S is a parameter needed for the boundary integral method. It should be a
% point(s) in a cavity/cavities away from the boundary/boundaries. 
S = 0 + 1i*0;

% Number of Panels for discretising the respective boundaries.
curve(1).nPanel = 64;
curve(2).nPanel = 44;

% Number of nodes; should be 16 or 32 (32 not implemented).
nNode = 16;
curve(1).nPoint = curve(1).nPanel * nNode;
curve(2).nPoint = curve(2).nPanel * nNode;

% Set plotFlag to 1 to plot results, otherwise 0.
plotFlag = 1;   

% The total number of boundaries, which we refer to as bodys in the co.
nBody = length(curve);



%% Discretise domain

% setupBoundary discretises the boundary with nNode Gaus-Legendre panels.
% It also produces the associated Gauss-Legendre weights, panel lengths.
% BIM_struct: see setup_boundary for explanation of elements.
% 
BIMstruct = setupBoundary(nNode, curve);
%
% Create a uniform grid for the Box B with M^2 points.
x = linspace(-LBox, LBox, M + 1);
x(end) = [];
[X, Y] = meshgrid(x,x);
xe = [X(:) Y(:)];

% Alapl is used for the boundary integral method for multiple connected domains.
% For identifying points as in or outside D we
% can set it to a zero vector of length nBody - 1. Not required if 
% nBody  = 1, i.e. no cavities.
BIMstruct.ALapl = 0;
% 
% idxbO is a boolean vector identifying points from xe as 
% in Omega (1) or outside (0).
 idxbO= inDomain(BIMstruct,curve,eps,xe);

% load('example2_800.mat')
% Points outside Omega, i.e. in E.
idxbE=~idxbO;

%% PDE setup
% Choose boundary conditions and right hand side.
[bc, f] = rhs.example2();
% Store general parameters and function handles in params.
params = struct('M',M,'bc',@(z) bc(z),...
    'f', @(z) f(z),...
    'plotFlag',plotFlag,'nBody',length(curve));




% Initialise PUX by precomputing all weights, interpolation matrix,
% distributing partitions centres and setting PUX parameters. All is stored
% in PUXstruct. Also an extension is given: fe. It has compact support, is
% equal to f in Omega. fe is a vector of length M^2.
[fe,PUXstruct] = setupPUX(M,xe,LBox,curve,params.f,BIMstruct.z,BIMstruct.arcL,BIMstruct.idxBody,PUXstruct,idxbO);


% Plot the extension and distribution of partitions.
if plotFlag
    plotExtension(PUXstruct,BIMstruct.z,BIMstruct.idxBody,xe,LBox,fe,params.f)
end





%% Particular solution
% Solve free-space with FFTs.

% uP is the particular solution on the box B with resolution M.
[uP,uHat, fHat] = freeSpacePoissonFourierCoeff(M,2*LBox,-fe);
uP = reshape(uP,[size(xe,1),1]);
% uHat are the Fourier coefficients for uP.
uHat = fftshift(uHat);


% Compute particular solution on the Gaus-Legendre points on all panels 
% along\Gamma.
uP_bdry = callNUFFT(uHat,M,BIMstruct.z,LBox,1e-14);


%% Homogeneuous solution

% Boundary Integral Equation, solve for density densLapl. 

[muLapl,ALapl] = dlpDensitySolver(sum(BIMstruct.nPanel)*BIMstruct.nNode,length(curve)-1,...
    BIMstruct.idxBody,BIMstruct.z,BIMstruct.dz, ...
    BIMstruct.d2z,BIMstruct.w,params.bc([real(BIMstruct.z) imag(BIMstruct.z)])-uP_bdry,S);

% Store density and A_Lapl (need name)
BIMstruct.mulapl = muLapl;
BIMstruct.ALapl = ALapl;

% Save evaluations points in complex form.
ze = [xe(idxbO,1)+1i*xe(idxbO,2)];


% Compute uH with standard quadrature.
uH = dlpDomainEval(muLapl, ALapl, S, ze,...
    BIMstruct.z, BIMstruct.dz, BIMstruct.w);

%% Apply special quadrature to treat near singularities.



for kBody = 1:nBody
idx = (1 + BIMstruct.idxBody(kBody)):BIMstruct.idxBody(kBody+1);
[uH] = dlpDomainEvalSpecQuad(uH, muLapl(idx), BIMstruct.nPanel(kBody), BIMstruct.z(idx), BIMstruct.dz(idx), BIMstruct.w(idx), ze, BIMstruct.zP((kBody + BIMstruct.idxBody(kBody)/BIMstruct.nNode):(BIMstruct.idxBody(kBody+1)/BIMstruct.nNode + kBody)),  IP1, IP2, W16, W32);
end



%% Error
% Create the analytical solution.
uSol = bc(xe);
% Allocate space for point wise error vector.
errAbsPntws = zeros(size(xe,1),1);
errAbsPntws(idxbO) = abs(uH+uP(idxbO)-uSol(idxbO));
errAbsPntws = errAbsPntws + eps;
%%
figure
surf(X,Y,reshape(log10(errAbsPntws),size(X)))
shading interp
colorbar
caxis([-16 1])
view(2)

maxError = max(errAbsPntws)/max(uSol(:));
relError = sqrt(sum(sum(errAbsPntws.^2))/sum(sum(abs(uSol).^2)));
disp(['Error on Grid = ', num2str(maxError)])

disp(['Relative error = ', num2str(relError)])
disp(' ')
errorL2 = relError;
errorMax = maxError;
end