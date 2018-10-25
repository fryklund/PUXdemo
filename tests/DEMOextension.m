%% DEMO: Extend function on simply/multiply connected domain.
% Notes:
% Boundary discretisation is mainly used to identify points as inside or
% outside Omega. 
%
close all;
clear all;
clc

%% Choose geometry and function to extend.

geometry_nr = 3;
function_nr = 1;

%% Construct Domain
% M: M^2 = Number of uniform grid points in Box

% LBox: Size of Box is [-LBox,LBox]^2. 
% n_Point: number of boundary points. In this demo mainly used to identify
% points in the interior/exterior of Omega.

switch geometry_nr
    case 1
           % One circle, with centre at (17/701, 5/439) and radius 1.
        LBox = 1.5;
        curve = diffCurve(curves.circle(17/701 + 1i*5/439,1));
        curve.nPoint = 1600;
    case 2
        % Two circles, with centre at (17/701, 5/439), with radius 1 and 0.7,
% respectively.
        
        LBox = 1.5;
        curve(1) = diffCurve(curves.circle(17/701 + 1i*5/439,1));
        curve(2) = diffCurve(curves.circleExterior(17/701 + 1i*5/439,0.7));
        curve(1).nPoint = 1600;
        curve(2).nPoint = 320;
        
    case 3
        % Two starfish/flower shaped domains.
        LBox = 1.8;
        curve(1) = diffCurve(curves.starfish(0.3,5));
        curve(2) = diffCurve(curves.starfishExterior(0.2,6));
        curve(1).nPoint = 1600;
        curve(2).nPoint = 320;

    otherwise 
        disp('Choose a geometry.')
        return;
end
    
        

nBody = length(curve);




%% Initialize domain

% Create uniform grid on Box.
M = 200;
xLin = linspace(-LBox, LBox, M+1);

% Make periodic
xLin(end) = [];
[X_Box, Y_Box] = meshgrid(xLin,xLin);
xe = [X_Box(:) Y_Box(:)];

% Discretize boundary and compute arc lengths.
z = []; % Set of boundary points
idxBody = zeros(nBody+1,1); % Vector pointers to z
arcL = zeros(nBody,1); % Vector of arc lengths.

for k = 1:nBody
    z = [z;curve(k).tau(linspace(0,2*pi,curve(k).nPoint)')];
    idxBody(k+1) = idxBody(k)+curve(k).nPoint;
    arcL(k) = integral(@(t) abs(curve(k).dtau(t)),0,2*pi);
end


%% List nodes that are inside or outside Omega.
% If point xe(i,:) is inside element then idxbO(i) = 1, oterhwise 0.

if nBody == 1
    idxbO = inpolygon(xe(:,1), xe(:,2), real(z), imag(z));
else
    inPolyBdry = [];
    for kBody = 1:nBody
        inPolyBdry = [inPolyBdry;NaN+1i*NaN;z((1+idxBody(kBody)):idxBody(kBody+1))];
    end
    inPolyBdry = [real(inPolyBdry(2:end)),imag(inPolyBdry(2:end))];
    idxbO = inpolygon(xe(:,1), xe(:,2), inPolyBdry(:,1), inPolyBdry(:,2));
end
% List containing the indices of points xe that are outside Omega.
idxbE =~idxbO;


%% Create function to extend
switch function_nr
    case 1
        % Gaussian
        [~,f] = rhs.DEMO_1();  
    case 2
        % Sinus
        [~,f] = rhs.example1();        
    otherwise 
        disp('Choose a function to extend.')
        return;
end

params = struct('M',M,'f',@(z) f(z),...
    'nBody',length(curve));
%% Set PUX parameters

% Width of interpolation basis, scales as 1/ep^2.
ep = 2;

% Sample data from uniform grid to create extension. 1 means use all
% available data, 2 means use every other point, 4 every fourth point etc.
% Use only even values.
coarsen = 1; 
% R_ratio sets the amount of overlap of the partitions. (1 = coverage, >1
% overlap). 2.5 is the recommended value.
R_ratio = 2.5;

% R_part is the partition radius.
R_p  = 0.4;

% regularity sets the C^k space the extension will belong to (k =
% regularity). If regularity = 0 it is set automatically.
regularity = 0; 

% Store values in struct.
PUXstruct.params = struct('ep',ep,'coarsen',coarsen,'R_ratio',R_ratio,'R_p',R_p,'regularity',regularity);


%% Init Pux
% Struct containing all information associated with the boundary.
% Uses z to check dist to bdry

[fe,PUXstruct] = setupPUX(M,xe,LBox,curve,params.f,z,arcL,idxBody,PUXstruct,idxbO);

plotExtension(PUXstruct,z,idxBody,xe,LBox,fe,params.f)
axis([-LBox LBox -LBox LBox])

figure(3)
surf(reshape(xe(:,1),[M,M]),reshape(xe(:,2),[M,M]),reshape(params.f(xe),[M,M]))
shading interp
view(2)
hold on
axis([-LBox LBox -LBox LBox])

for kBody = 1:length(curve)
    idx = (1+idxBody(kBody)):idxBody(kBody+1);
    plot3(real(z(idx)),imag(z(idx)),params.f([real(z(idx)),imag(z(idx))]),'k-','LineWidth',2)
end
fig3_cb = colorbar;

figure(2)
caxis(fig3_cb.Limits)
axis equal
axis([-LBox LBox -LBox LBox])

