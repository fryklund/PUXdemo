function [] = plotExtension(PUXstruct,z,idxBody,xe,LBox,fe,f)
nBody = length(PUXstruct.curve);
M = sqrt(size(xe,1));
figure(1)
% pcolor(reshape(Domain_struct.xe(:,1),[M,M]),reshape(Domain_struct.xe(:,2),[M,M]),reshape(fe,[M,M]))
% shading interp


t = linspace(0,2*pi,100);
hold on
axis equal
for kBody = 1:nBody
    idx = (1+idxBody(kBody)):idxBody(kBody+1);
   plot(real(z(idx)),imag(z(idx)),'k-','LineWidth',2)
    plot3(xe(PUXstruct.curve(kBody).idxInterpC,1),xe(PUXstruct.curve(kBody).idxInterpC,2),0*xe(PUXstruct.curve(kBody).idxInterpC,1),'k*')
    for i = 1:PUXstruct.curve(kBody).nInterp
        part = xe(PUXstruct.curve(kBody).idxInterpC(i),1)+1i*xe(PUXstruct.curve(kBody).idxInterpC(i),2) + exp(1i*t)*PUXstruct.curve(kBody).R_Interp;
        plot3(real(part),imag(part),0*part,'k-')
    end      
    
end
%%
for kBody = 1:nBody
    plot3(xe(PUXstruct.curve(kBody).idxZeroC,1),xe(PUXstruct.curve(kBody).idxZeroC,2),0*xe(PUXstruct.curve(kBody).idxZeroC,1),'bo')
    for i = 1:PUXstruct.curve(kBody).nZero   
        part = xe(PUXstruct.curve(kBody).idxZeroC(i),1) + 1i*xe(PUXstruct.curve(kBody).idxZeroC(i),2)+ exp(1i*t)*PUXstruct.curve(kBody).R_Zero(i);
        plot3(real(part),imag(part),0*part,'b--')
    end
    
end
hold off
axis([-LBox LBox -LBox LBox])
figure(2)
surf(reshape(xe(:,1),[M,M]),reshape(xe(:,2),[M,M]),reshape(fe,[M,M]))
shading interp
view(2)
hold on


for kBody = 1:nBody
    idx = (1+idxBody(kBody)):idxBody(kBody+1);
    plot3(real(z(idx)),imag(z(idx)),f([real(z(idx)),imag(z(idx))]),'k-','LineWidth',2)
end
colorbar


