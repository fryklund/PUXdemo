function [idxbO] = inDomain(BIMstruct,curve,S,xe)



ze = xe(:,1) + 1i*xe(:,2);
inner = dlpDomainEval(ones(size(BIMstruct.z)), BIMstruct.ALapl, S, ze,...
    BIMstruct.z, BIMstruct.dz, BIMstruct.w);



load 'IP1632.mat'
load 'glW.mat' %read in GL 16 and 32 weights

for kBody = 1:length(curve)
    idx = (1 + BIMstruct.idxBody(kBody)):BIMstruct.idxBody(kBody+1);
    [inner] = dlpDomainEvalSpecQuad(inner, ones(length(idx),1), BIMstruct.nPanel(kBody), BIMstruct.z(idx), BIMstruct.dz(idx), BIMstruct.w(idx), ze, BIMstruct.zP((kBody + BIMstruct.idxBody(kBody)/BIMstruct.nNode):(BIMstruct.idxBody(kBody+1)/BIMstruct.nNode + kBody)),  IP1, IP2, W16, W32);
end
idxbO = false(size(inner));

idxbO(abs(inner-1)< 1e-12) = true;
end