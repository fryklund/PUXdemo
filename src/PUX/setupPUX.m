function [fe,PUXstruct] ...
    = setupPUX(M,xe,LBox,curve,f,z,arcL,idxBody,PUXstruct,idxbO)
%% ====================     PUX     =======================
% INPUT:
%   M:
%       M^2 number of points on uniform grid.
%   xe:
%       [M^2,2] matrix containing the uniform grid     
%   LBox:
%       Size of box B = [-LBox,LBox]^2
%   curve:
%       A structure containing functions corresponding to the
%       parametrization of the curve. All fields have the input t in [0,2pi],
%       The structure contains the following fields:
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
%   f:
%       Function to extend. Vector of function values at points in xe.
%   z: 
%       Complex. Discretization of boundary/boundaries. 
%   arcL:
%       Arc length of boundary/boundaries. 
%   idxBody:
%       A vector containing indices used to extract boundary specific
%       elements.
%   PUXstruct:
%       A structure containing parameters for PUX. Will store PUX data.
%       As input it contains only the field params, which contains the
%       following fields.
%       ep: 
%           Shape parameter for width Gaussian interpolation basis.
%           Scales as 1/ep. 
%       coarsen:
%           Sets how to sample the evaluation grid to obtain the
%           interpolation grid. Set to 1 if they should be equal.    
%       R_ratio:
%           Amount of partition overlap (1 = coverage, >1 overlap)
%       R_p:
%           Partition radius.
%       Regularity:
%           Integer that sets the regularity of weight functions. Regularity 1
%           to 5 are implemented. If 0 then regularity is set automatically to
%           for near optimal choice.
%   idxbO: OPTIONAL, indices for points in xe thar are in Omega.
% OUTPUT:
%   fe:
%       Extension of f at grid points in xe.
%   PUXstruct:
%       Struct passed as input has been updated with the following fields,
%       with the corresponding element.
%       PUXstruct.data
%           A: The mapping of data from nonuniform points to uniform, see
%           PUX-paper.
%           W: Matrix with partition of unity weights.
%           If: Indices for which points on uniform grid that are the partitions
%           extension points.
%           Jf: Indices for which points on uniform grid belong to a partition.
%       PUXstruct.curve(k), k =1:nBody
%           R_part: Radius of interpolation partitions (scalar).
%           R_Zero: Radius of zero partitions (vector).
%           nInterp: Number of interpolation partitions.
%           nZero: Number of zero partitions.
%           idxInterpC: indices for interpolation partions' centres in xe.
%           idxZeroC: indices for zero partions' centres in xe.


% LOCAL VARIABLES:
%   nBody: Number of boundaries.
%   R_ZeroThreshold: If a zero partition's radius is less than
%   R_ZeroThreshold* R_part then it is removed.
%   distR = Distance between partition centers.
%   nRBFC =  Number of Vogel points, i.e. RBF-centres.
%   P =  Measure of point density
%   nBoxPnts = Total number of points in uniform grid. 
%   h = Grid resolution
%   nPart = The total number of partitions. 

%   Comment on naming: We use the following abbrevations: idx =
%   index/indices, idxb = index/indices as boolean (1 for being in a set,
%   otherwise 0), Interp = interpolation partition, Zero = zero partition,
%   P = partition, C = centre (used as InterpC for interpolation partitions' centres),
%  n* = number of *, dist = distance, O = Omega (the domain/domains where
%   f is known, E = Extension domain (i.e. complement to Omega with respect
%   to surrounding box), col = collocation, Val = values from evaluation.
%   Unless specified otherwise, idx in combination with one set, e.g.
%   idxInterp means the intersection of uniform points in Interp (i.e.
%   interpolation partitions) and the uniform grid, i.e. points in xe. For
%   example idxInterpE means points in xe that are in E and the
%   interpolation partitions. A single partition is specified as Interp_j. Local indices are denoted
%   idxLoc.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set and extreact parameters.
nBody = length(curve);
R_ZeroThreshold = 0.8;
coarsen = PUXstruct.params.coarsen; 
R_ratio = PUXstruct.params.R_ratio; 
R_part=PUXstruct.params.R_p;
distR = R_part/ R_ratio;
P = R_part * M / (2 * LBox);
nRBFC = round(min(0.8*pi/4*P.^2,4*P));
nBoxPnts = size(xe,1);
nInterp = (ceil(arcL./(2*distR))+1);
h = abs(xe(2,2)-xe(1,2));
nPart = sum(nInterp); % nZero are added later.

%% Print some parameters
disp(['P = ' num2str(P)]);
disp(['Number of Vogel nodes:' num2str(nRBFC)])
disp(['Rp = ' num2str(R_part)]);


%% Create coarser grid for interpolation
if (coarsen > 1)
    [I,J] = meshgrid(1:coarsen:M);
    i_crs = I+M*(J-1);
    i_crs = i_crs(:);
end
%% List nodes that are inside or outside Omega.
% The boundary of Omega is given by complex z.

% Mark the points in xe as inside or outside Omega. If point xe(i,:) is
% inside, then idxbO(i) = 1, oterhwise 0.
if nargin <10
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
end

% List containing the indices of points xe that are outside Omega.
idxbE=~idxbO;

%% Find values in collocation points
fcol = zeros(nBoxPnts,1);
fcol(idxbO) = f(xe(idxbO,:));

%% Find centers for interpolation partitions
% We find the centres for the nInterp interpolation partitions InterpPart by
% sampling uniformly with respect to arclength along the boundary. For each
% such point the closest uniform grid point inside Omega is found and set
% to be a partition centre. Since 2*nInterp zero partitions centres are used we
% can sample them at the same time.


% Create a list of centres for interpolation partitions. 
InterpCList = zeros(sum(nInterp),2);
InterpCListCounter = 0;
% Loop over boundaries to distribute partitions.
for kBody = 1:nBody
    % Create new uniformly distributed points w.r.t. arclength.
    [ZeroCnonuniform,~,~] = Traparcldisc(curve(kBody),...
        nInterp(kBody)*2);
    InterpCnonuniform = ...
        [real(ZeroCnonuniform(1:2:end)),...
        imag(ZeroCnonuniform(1:2:end))];
    InterpCnonuniform = (InterpCnonuniform + LBox)/h;
    idxInterpC = zeros(nInterp(kBody),1);
    for j = 1:nInterp(kBody)
        floornonuniform = floor(InterpCnonuniform(j,:));
        ceilNonUniform = ceil(InterpCnonuniform(j,:));
        idxInterpCuniform = [floornonuniform(1)*M+...
            floornonuniform(2)+1;...
            floornonuniform(1)*M+...
            ceilNonUniform(2)+1;...
            ceilNonUniform(1)*M+...
            floornonuniform(2)+1;...
            ceilNonUniform(1)*M+...
            ceilNonUniform(2)+1];
        % Requires that the box E is centered at the origin, without taking
        % periocidity into acccount.
        idxInterpCcandidates = find(idxbO(idxInterpCuniform)); 
        idxInterpC(j) = idxInterpCuniform(idxInterpCcandidates(1));
    end
    
    R_Interp = R_part;
    
    InterpCList((1+InterpCListCounter):(InterpCListCounter+length(idxInterpC)),:) = xe(idxInterpC,:);
    InterpCListCounter = InterpCListCounter + length(idxInterpC);
    %% Create zero partitions
    
    nZero = nInterp(kBody)*2;
    
    
    %dZeroCnonuniform contains the values of the derivative of the
    %paramterization of the curve evaluated at the nZero points.
    [dZeroCnonuniform,~] = fftDiff(ZeroCnonuniform ,[1,nZero]);
    ZeroCnonuniform = [real(ZeroCnonuniform),imag(ZeroCnonuniform)];
    if kBody == 1
        ZeroCnonuniformNormal = -1i*dZeroCnonuniform./abs(dZeroCnonuniform);
    else
        ZeroCnonuniformNormal = -1i*dZeroCnonuniform./abs(dZeroCnonuniform);
    end
    ZeroCnonuniform = [ZeroCnonuniform(:,1) + (R_Interp + h/2)*real(ZeroCnonuniformNormal), ZeroCnonuniform(:,2) + (R_Interp + h/2)*imag(ZeroCnonuniformNormal)];
    
    
    idxZeroC = round((ZeroCnonuniform+LBox)/h);
    idxZeroC = idxZeroC(:,1) * M + idxZeroC(:,2) + 1;
    
    % Some partition centres might end inside Omega, happens in concave ares. These are removed.
    idxZeroC = idxZeroC(logical(idxbE(idxZeroC)));
    
    
    % Set radius for zero partitions, remove those that are too small.
    [~,distBdry2ZeroC] = knnsearch([real(z) imag(z)],xe(idxZeroC,:));
    % [~,Dist_bdry2ZeroPartCentres] = knnsearch([real(z(:,1)) imag(z(:,1))],xe(idxZeroC,:));
    idxZeroKeep = find(distBdry2ZeroC>R_Interp*R_ZeroThreshold);
    R_Zero = distBdry2ZeroC(idxZeroKeep);
    idxZeroC = idxZeroC(idxZeroKeep);
    % If a zero-partition has greater than R_Ip radius, set its radius to R_Ip.
    R_Zero(logical(R_Zero>R_Interp)) = R_Interp;
    nZero = size(idxZeroC,1);
    nPart = nPart + nZero;
    
    %% Not included: 'idx_xeOutOmegaInInterpPart',idx_xeOutOmegaInInterpPart
    
    PUXstruct.curve(kBody) = struct('R_Interp',R_Interp,'R_Zero',R_Zero,'nInterp',nInterp(kBody),'nZero',nZero,...
        'idxInterpC',idxInterpC,'idxZeroC',idxZeroC);
    % Create Vogel-points. These will be used as RBF-centers.
    
    
end
%% Create a NeighborSearcher object for K-nearest neighbors search.
NeighborSearcher = createns(xe);

%% Create interpolation stencil

% Locate index in xe for points withing R_Interp of InterpPart.
% Sufficient to do this for one interpolation partition and then shift.
[idxInterpcell,distGrid2InterpC] = rangesearch(NeighborSearcher,...
    xe(PUXstruct.curve(1).idxInterpC(1),:),...
    PUXstruct.curve(1).R_Interp);
[idxInterpStencil, idxSortAscen] = sort(idxInterpcell{1}-PUXstruct.curve(1).idxInterpC(1));
distGrid2InterpCstencil = distGrid2InterpC{1};
distGrid2InterpCstencil = distGrid2InterpCstencil(idxSortAscen)';
nPntsInterp = length(idxInterpStencil);

%% Choose RBF to use as building block for weight function.
if PUXstruct.params.regularity == 0
    PUXstruct.params.regularity= min(floor(sqrt(P)-0.9),5); % test
end
switch PUXstruct.params.regularity
    case 1
        disp('Shepard with Wu1')
        Phi = @(r) 1/2*(1-r).^2.*(2+r); % Wu C0(R3)
    case 2
        disp('Shepard with Wu2')
        Phi = @(r) 1/8*(1-r).^3.*(8+9*r+3*r.^2); % Wu C0(R5)
    case 3
        disp('Shepard with Wu3')
        Phi = @(r) 1/4*(1-r).^4.*(4+16*r+12*r.^2+3*r.^3); % Wu C2(R3)
    case 4
        disp('Shepard with Wu4')
        Phi = @(r) 1/8*(1-r).^5.*(8+40*r+48*r.^2+25*r.^3+5*r.^4); % Wu C2(R5)
    case 5
        disp('Shepard with Wu5')
        Phi = @(r) 1/6*(1-r).^6.*(6+36*r+82*r.^2+72*r.^3+30*r.^4+5*r.^5); % Wu C4(R3)
    otherwise
        disp('Shepard with Wendland C2')
        Phi = @(r) (1-r).^4.*(4*r+1); % Wendland C2(R3)
end

%% Distribute RBF-centres
% Vogel point distribution.
RBFC = PUXstruct.curve(1).R_Interp *...
    repmat(sqrt(1:nRBFC(1))'/sqrt(nRBFC(1)),1,2).*...
    [cos((1:nRBFC(1))'*pi*(3-sqrt(5))),...
    sin((1:nRBFC(1))'*pi*(3-sqrt(5)))];
%% Evaluate an RBF associated with an interpolation partition once, then reuse.
RBFInterpStencilVal = Phi(distGrid2InterpCstencil/R_Interp);
%% Set up for extension
% Total number of weight function evaluations.
nWeightVal = sum(nInterp) * nPntsInterp * 3;
% Crude approximation of number of points we extend to. Used to allocate
% space.
nExtPnts = round(sum(nInterp) * nPntsInterp*2/3);

% Vectors for data used for building extension. The names I, J and S are
% based on the input for MATLAB's sparse:
% S = sparse(i,j,s,m,n) uses vectors i, j, and s to generate an
%     m-by-n sparse matrix such that S(i(k),j(k)) = s(k)
If = zeros(nExtPnts,1); Jf = zeros(nExtPnts,1); Sf = zeros(nExtPnts,1);
IW = zeros(nWeightVal,1); JW = zeros(nWeightVal,1); SW = zeros(nWeightVal,1);
% Counters used to fill in If, Jf, Sf, IW, JW and SW.
idxEndf = 0;
idxEndW = 0;
%% RBF-Direct
% Uncomment to use RBF-Direct.
% xe1 = xe(idx_inIp + idx_Ip(1),:);
%     xv1 = x_rbf+repmat(Ip(1,:),n_v,1);
%
%     [x1,x2] = ndgrid(xe1(:,1),xv1(:,1));
%     [y1,y2] = ndgrid(xe1(:,2),xv1(:,2));
%     r2 = (x1-x2).^2 + (y1-y2).^2;
%     A = exp(-params.ep^2*r2); xe(idx_xeInInterpPart_stencil + PUX_struct.curve(1).idxInterpC(1),:)
%     PUXstruct.general.A = A;
%% RBF--QR
% Create PUX-matrix for the stencil interpolation partition, using RBF-QR.
A = RBF_QR_diffmat_2D('1',...
    xe(idxInterpStencil + PUXstruct.curve(1).idxInterpC(1),:),...
    [RBFC(:,1)+xe(PUXstruct.curve(1).idxInterpC(1),1),...
    RBFC(:,2)+xe(PUXstruct.curve(1).idxInterpC(1),2)],...
    PUXstruct.params.ep, 2);
%% Loop over interpolation partitions
msgid = 'MATLAB:nearlySingularMatrix';
warning('off',msgid);
msgid2 = 'MATLAB:rankDeficientMatrix';
warning('off',msgid2);
% minData stores the number of interpolation points for the
% interpolation partitions with the fewest. Used to determine if the least
% squares systems are sufficiently "overdetermined".
minData = 10^16;
% Same as minData, but stores the maximum amount of interpolation
% points. Currently not used.
maxData = 0;
% Need to keep track of which index to start at for each partition.
idxJfstart = 0;
for kBody = 1:nBody
    for j = 1:PUXstruct.curve(kBody).nInterp
        % Global indices
        idxInterp_j = PUXstruct.curve(kBody).idxInterpC(j) + idxInterpStencil; % Indices of point in xe within R_Ip of parition center Ip(i).
        idxInterp_jO = idxInterp_j(idxbO(idxInterp_j));
        idxInterp_jE = idxInterp_j(~idxbO(idxInterp_j));
        nLocalV = length(idxInterp_jE);       
        % Local indices
        idxLocInterp_jO = find(idxbO(idxInterp_j));
        idxLocInterp_jE = find(~idxbO(idxInterp_j));
        if (coarsen > 1) %TODO: Do once and reuse.
            % Coarser grid. Samples points set by coarsen.
            [ic,ij] = intersect(idxInterp_j(idxbO(idxInterp_j)),i_crs); % global
            idxInterp_jO = ic(idxbO(ic));  % global
            idxLocInterp_jO = idxLocInterp_jO(ij);
        end
        % Solve overdetermined system to obtain f in RBF centres.
        fRBFC = A(idxLocInterp_jO,:)\fcol(idxInterp_jO);
        maxData = max(length(idxInterp_jO),maxData);
        minData = min(length(idxInterp_jO),minData);
        If(idxEndf+1:idxEndf+nLocalV) = idxInterp_jE;
        Jf(idxEndf+1:idxEndf+nLocalV) = j + idxJfstart;
        Sf(idxEndf+1:idxEndf+nLocalV) = A(idxLocInterp_jE,:) * fRBFC;
        
        % Uncomment to use analytic values for extension.
        %   Sf(idx_endf+1:idx_endf+n_local_evals) =  params.f(xe(idx_eval_ip,:)); %    Use actual f for extension.
        % warning(['No extrapolation. Analytical f evaluated in E.']);
        
        % Save contribution to weight function
        IW(idxEndW+1:idxEndW+nPntsInterp) = idxInterp_j;
        JW(idxEndW+1:idxEndW+nPntsInterp) = j + idxJfstart;
        SW(idxEndW+1:idxEndW+nPntsInterp) = RBFInterpStencilVal;
        idxEndf = idxEndf+nLocalV;
        idxEndW = idxEndW+nPntsInterp;       
    end
    idxJfstart = idxJfstart + PUXstruct.curve(kBody).nInterp;
end
If = If(1:idxEndf); Jf = Jf(1:idxEndf); Sf = Sf(1:idxEndf);
idxInterpE = sort(unique(If));
Imat = sparse(If,Jf,Sf,nBoxPnts,sum(nPart));
warning('on',msgid)
warning('on',msgid2)
beta_min = minData/nRBFC;
if beta_min<2
    warning(['Least-squares oversampling less than two (' num2str(beta_min) ')!']);
end
disp(['beta_min: ' num2str(beta_min)])
%% Calculate contribution to weight function for zero partitions
idxJWstart_idx = sum(nInterp);
for kBody = 1:nBody
    for j = 1:PUXstruct.curve(kBody).nZero
        idxZeroPart_j = PUXstruct.curve(kBody).idxZeroC(j) + idxInterpStencil;
        % Some partitions fall outside box [-LBox,LBox]^2. Prune uniform
        % grid accordingly.
        if  max(idxZeroPart_j)> nBoxPnts || min(idxZeroPart_j)< 0
            idx_insideB = intersect(find(idxZeroPart_j < nBoxPnts+1),find(idxZeroPart_j > 0));
            nPntsZero_j = length(idx_insideB);
            distGird2ZeroC_j = distGrid2InterpCstencil(idx_insideB);
            IW(idxEndW+1:idxEndW+nPntsZero_j) = idxZeroPart_j(idx_insideB);
            JW(idxEndW+1:idxEndW+nPntsZero_j) = j+idxJWstart_idx;
            SW(idxEndW+1:idxEndW+nPntsZero_j) = Phi(distGird2ZeroC_j/PUXstruct.curve(kBody).R_Zero(j));
            idxEndW = idxEndW+nPntsZero_j;
            % If the radius of a zero partition is not R_interp.
        elseif norm(PUXstruct.curve(kBody).R_Zero(j)- PUXstruct.curve(kBody).R_Interp) > 1e-14
            distGird2ZeroC_j = find(distGrid2InterpCstencil<=PUXstruct.curve(kBody).R_Zero(j));
            idxZeroPart_j = idxZeroPart_j(distGird2ZeroC_j);
            distGird2ZeroC_j = distGrid2InterpCstencil(distGird2ZeroC_j);
            nPntsZero_j = length(distGird2ZeroC_j);
            IW(idxEndW+1:idxEndW+nPntsZero_j) = idxZeroPart_j;
            JW(idxEndW+1:idxEndW+nPntsZero_j) = j+idxJWstart_idx;
            SW(idxEndW+1:idxEndW+nPntsZero_j) = Phi(distGird2ZeroC_j/PUXstruct.curve(kBody).R_Zero(j));
            idxEndW = idxEndW+nPntsZero_j;
        else
            IW(idxEndW+1:idxEndW+nPntsInterp) = PUXstruct.curve(kBody).idxZeroC(j) + idxInterpStencil;
            JW(idxEndW+1:idxEndW+nPntsInterp) = j+idxJWstart_idx;
            SW(idxEndW+1:idxEndW+nPntsInterp) = RBFInterpStencilVal;
            idxEndW = idxEndW+nPntsInterp;
        end
        
    end
    % Update index.
    idxJfstart = idxJfstart + PUXstruct.curve(kBody).nZero;
end
%% Build weight matrix W
% Remove zero-valued elements, these occur because all paritions do not
% hold nPntsInterpPart points.
IW = IW(1:idxEndW); JW = JW(1:idxEndW); SW = SW(1:idxEndW);
W = sparse(IW,JW,SW,nBoxPnts,nPart);
W = W(idxInterpE,:);
Imat = Imat(idxInterpE,:);
sumW = sum(W,2);
% W contains the Partition of Unity weights. Can be resused.
W = spdiags(1./sumW,0,length(idxInterpE),length(idxInterpE))*W;
%% Combine local extensions into global extension
fe = zeros(nBoxPnts,1);
% Combine local extensions with partition of unity weights.
fe(idxInterpE) = full(sum(Imat.*W,2));
% Enforce fe = f on Omega.
fe(idxbO) = fcol(idxbO);
PUXstruct.data = struct('A',A,'W',W,'If',If,'Jf',Jf,'idxInterpStencil',idxInterpStencil,'idxbO',idxbO);
end