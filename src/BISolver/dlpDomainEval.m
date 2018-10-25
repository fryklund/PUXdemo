function uH = dlpDomainEval(muLapl,ALapl, S, ze, z, dz, w)

uH = zeros(size(ze));
% for j=1:length(z)
%     for k=1:N
%         if z(j) ~= zDrops(k)
%              u(j) = u(j) + 1/(2*pi)*wDrops(k)*mu_lapl(k)*imag(zpDrops(k)/(zDrops(k)-z(j)));
%         else
%             u(j) = u(j) + 1/(2*pi)*wDrops(k)*mu_lapl(k)*imag(zppDrops(k)/(zpDrops(k)));
%         end
%     end
% end

parfor j=1:length(ze)
   uH(j) = sum(muLapl.*w.*imag(dz./(z-ze(j)))) + 2*pi*sum(log(abs(ze(j)-S)).*ALapl); 
end
uH = 1/(2*pi)*uH;

end