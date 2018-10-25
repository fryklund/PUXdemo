function [fe] = compute_function_extension(xe,f,PUX_struct)
%% Extraxt data from struct

A_PUX = PUX_struct.data.A_PUX;
W_PUX = PUX_struct.data.W_PUX;
If = PUX_struct.data.If;
Jf = PUX_struct.data.Jf;
n_Boxpnts = size(xe,1);
nBody = length(PUX_struct.curve);
n_Part = 0;


idx_xeOutOmegaInInterpPart = sort(unique(PUX_struct.data.If));
idx_xeInOmega_boolean = PUX_struct.data.idx_xeInOmega_boolean;
idx_xeInInterpPart_stencil = PUX_struct.data.idx_xeInInterpPart_stencil;
coarsen = PUX_struct.params.coarsen;

%% Find values in collocation points
f_collocation = zeros(n_Boxpnts,1);
% f_collocation(idx_xeInOmega_boolean) = f(xe(idx_xeInOmega_boolean,:));
f_collocation(idx_xeInOmega_boolean) = f(idx_xeInOmega_boolean);

%% 
% Create coarser grid for interpolation
if (coarsen > 1)
    M = sqrt(size(xe,1));
    [I,J] = meshgrid(1:coarsen:M);
    i_crs = I+M*(J-1);
    i_crs = i_crs(:);
end


% Vector for data used for building extension. The names S is
% based on the input for MATLABS sparse:
% S = sparse(i,j,s,m,n) uses vectors i, j, and s to generate an
%     m-by-n sparse matrix such that S(i(k),j(k)) = s(k)

% Sf stores information about the local extensions.
Sf = zeros(size(Jf));


% Counters used to fill in Sf.
idx_endf = 0;
msgid = 'MATLAB:nearlySingularMatrix';
warning('off',msgid);
msgid2 = 'MATLAB:rankDeficientMatrix';
warning('off',msgid2);

for kBody = 1:nBody
for i = 1:PUX_struct.curve(kBody).n_InterpPart
    % Global indices
    idx_xeInInterpPart_i = PUX_struct.curve(kBody).idx_InterpPartCentres(i) + idx_xeInInterpPart_stencil; % Indices of point in xe within R_Ip of parition center Ip(i).
    idx_xeInOmegaInInterpPart_i = idx_xeInInterpPart_i(idx_xeInOmega_boolean(idx_xeInInterpPart_i));
    idx_xeInInterpPartOutOmega = idx_xeInInterpPart_i(~idx_xeInOmega_boolean(idx_xeInInterpPart_i));
    n_local_evals = length(idx_xeInInterpPartOutOmega);
    % Local indices
    idx_local_xeInOmega = find(idx_xeInOmega_boolean(idx_xeInInterpPart_i));
    idx_local_xeOutOmega = find(~idx_xeInOmega_boolean(idx_xeInInterpPart_i));

    if (coarsen > 1)
        % Coarser grid. Samples points set by coarsen.
        [ic,ij] = intersect(idx_xeInInterpPart_i(idx_xeInOmega_boolean(idx_xeInInterpPart_i)),i_crs); % global
        idx_xeInOmegaInInterpPart_i = ic(idx_xeInOmega_boolean(ic));  % global
        idx_local_xeInOmega = idx_local_xeInOmega(ij);
    end
    % Solve overdetermined system to obtain f in RBF centres.
    f_RBFCentres = A_PUX(idx_local_xeInOmega,:)\f_collocation(idx_xeInOmegaInInterpPart_i);
    Sf(idx_endf+1:idx_endf+n_local_evals) = A_PUX(idx_local_xeOutOmega,:) * f_RBFCentres; %#ok<FNDSB> %
% 
%   Sf(idx_endf+1:idx_endf+n_local_evals) =  f(idx_xeInInterpPartOutOmega,:); %    Use actual f for extension.
% warning(['No extrapolation. Analytical f evaluated in E.']);
% 
%     
    idx_endf = idx_endf+n_local_evals;
    
end
n_Part = n_Part + PUX_struct.curve(kBody).n_InterpPart+PUX_struct.curve(kBody).n_ZeroPart;
end
% cell_idx_xeOutOmegaInInterpPart is obtained by range search from 
%partition centers on uniform grid. Thus they are not the same amount as 
%idx_local_xeOutOmega_boolean. idx_local_xeOutOmega_boolean
% comes from inpolygon. 
Sf = Sf(1:idx_endf);
Imat = sparse(If,Jf,Sf,n_Boxpnts,n_Part);


Imat = Imat(idx_xeOutOmegaInInterpPart,:);


%% Combine local extensions into global extension
fe = zeros(n_Boxpnts,1);

fe(idx_xeOutOmegaInInterpPart) = full(sum(Imat.*W_PUX,2));


fe(idx_xeInOmega_boolean) = f_collocation(idx_xeInOmega_boolean);

end







