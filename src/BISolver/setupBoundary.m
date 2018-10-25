function [BIMstruct]= setupBoundary(nNode, curve)
% INPUTS:
%   nNode:
%       Number of Gauss-Legendre nodes on a panel.
%   curves:
%       a structure array containing functions corresponding to the
%       parametrization of the curve. All fields have the input t, which is assumed to go from 0 to 2*pi. The structure contains the following
%       fields:
%
%       curves.tau:
%           The parametrization of the curve
%       curves.dtau:
%           The derivative of the parametrization
%       curves.d2tau:
%           The second derivative of the parametrization
%       curves.normal:
%           The inward normal at point curves.tau(t) to the curve given by tau.
%       curve.tangent:
%           The tangent at point curves.tau(t) to the curve given by tau.
%       curve.kappa: (Not yet implemented)
%           The curvature at point curves.tau(t) to the curve given by tau.
%       curve.nPanel:
%           Number of panels.
% OUTPUTS: 
% BIMstruct, which contains the following elements:
%
%   t:  parameter values at Gauss-Legendre nodes
%       
%   z, dz, d2z:
%       a nPoints x nBody matrix consisting of the points along the curve,
%       the first and second derivatives with respect to the
%       parametrization, respectively. Boundary points are equispaced
%       in the parameter, and listed in a counter clockwise
%       direction. The orientation of dz is determined by the relation
%       of the curve with respect to the domain.
%   zP:
%       Endpoints of panels on curve
%   ds:
%       |dz/dt|
%   Nz:
%       Outward pointing normal
%   kappa:
%       Curvature
%   T:
%       canonical Gauss-Legendre nodes
%   W:
%       canonical Gauss-Legendre weights
%   idxBody:
%       vector of indices used to access the elements from the Gauss--Legendre points and panels. Length nBody + 1.
%


%% =====================     BIM setup       ============================
% Create panels and Gauss-Legendre points for each panels.

nBody = length(curve);
idxBody = zeros(nBody+1,1);
nPanel = zeros(nBody,1);

if nBody > 1
    nPoint = 0;
for kBody = 2:nBody+1
    idxBody(kBody) = curve(kBody-1).nPoint+idxBody(kBody-1);
    nPoint = nPoint + curve(kBody-1).nPoint;
    nPanel(kBody-1) = curve(kBody-1).nPanel;
end
else
    nPoint = curve(1).nPanel * nNode;
    idxBody(1) = 0;
    idxBody(2) = nPoint;
    nPanel = curve(1).nPanel;
end


t = zeros(nPoint,1);
[T, W] = GaussLegendre(nNode);
w = zeros(nPoint,1);
z = zeros(nPoint,1);
dz = zeros(nPoint,1);
d2z = zeros(nPoint,1);
zP = zeros(nPoint/nNode+nBody, 1);
Nz = zeros(nPoint,1);

arcLNew = zeros(nBody,1);
for kBody = 1:nBody
    dt = 2*pi/curve(kBody).nPanel;
    w_kBody = 0.5*dt*W;   % Regular quadrature weights
    for iPanel = 1:curve(kBody).nPanel
        idx = (1:nNode) + (iPanel-1)*nNode + idxBody(kBody);
        ta = (iPanel-1)*dt;
        t(idx) = (ta + 0.5*dt*(T + 1));
        z(idx) = curve(kBody).tau(t(idx));
        dz(idx)= curve(kBody).dtau(t(idx));
        d2z(idx)= curve(kBody).d2tau(t(idx));
        zP(iPanel + idxBody(kBody)/nNode+(kBody-1)) = curve(kBody).tau(ta);
        w(idx) = w_kBody;
    end
    arcLNew(kBody) = sum(w((1+idxBody(kBody)):idxBody(kBody+1)).* abs(abs(dz((1+idxBody(kBody)):idxBody(kBody+1)))));
    zP(iPanel + 1 + idxBody(kBody)/nNode+(kBody-1)) = zP(1 + (kBody-1) + idxBody(kBody)/nNode);    
end
ds = abs(dz);
Nz(:)= -1i*dz./abs(dz);% outward normal
kappa = -imag(dz.*conj(d2z))./ds.^3;    % curvature


LGammaP = zeros(nPoint/nNode, 1);
arcL_disc = zeros(nBody,1);
arcL = zeros(nBody,1);
[T_1000,W_1000] = GLinterval(0,2*pi,1000);
for kBody = 1:nBody
    for iPanel = 1: curve(kBody).nPanel
        i = (iPanel-1)*nNode + (1: nNode);
        LGammaP(iPanel+idxBody(kBody)/nNode) = sum(ds(i+idxBody(kBody)).*w((1:nNode) + nNode*(kBody-1)));
    end
    arcL_disc(kBody) = sum(LGammaP((1+idxBody(kBody)/nNode):idxBody(kBody+1)/nNode));
    arcL(kBody) = W_1000'*abs(curve(kBody).dtau(T_1000));
end









BIMstruct = struct('z',z,'dz',dz,'d2z',d2z,'zP',zP,'ds',ds,'Nz',Nz,'kappa',kappa,'t',t,'w',w,'W',W,'T',T,'LGammaP',LGammaP,'arcL_disc',arcL_disc,'arcL',arcL,'nNode',nNode,'idxBody',idxBody,'nPanel',nPanel);
end







