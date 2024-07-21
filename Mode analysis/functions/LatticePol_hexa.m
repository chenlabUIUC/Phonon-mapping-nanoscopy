function [Lp,Lm]=LatticePolAng(x,l0)
%x=[k]
%Secondary Constants
x=abs(x);
k=x(1);
ang = 60; % deg
%Vectors
a1=[1,0];
a2=[cosd(ang),sind(ang)];
a3=[cosd(ang*2),sind(ang*2)];
A1=a1*l0; A2=a2*l0; A3=a3*l0;

%Matrix
M=@(kx,ky,a,b) 2*k*(a1(a)*a1(b).*(1-cos(kx*A1(1)+ky*A1(2)))+...
                    a2(a)*a2(b).*(1-cos(kx*A2(1)+ky*A2(2)))+...
                    a3(a)*a3(b).*(1-cos(kx*A3(1)+ky*A3(2))));

SQ= @(k1,k2) sqrt(((M(k1,k2,1,1)+M(k1,k2,2,2))./2).^2-(M(k1,k2,1,1).*M(k1,k2,2,2)-M(k1,k2,1,2).*M(k1,k2,2,1)));

yp= @(k1,k2) (((M(k1,k2,2,2)-M(k1,k2,1,1))./2)+SQ(k1,k2))./sqrt(M(k1,k2,1,2).*M(k1,k2,2,1)+(((M(k1,k2,2,2)-M(k1,k2,1,1))./2)+SQ(k1,k2)).^2);
ym= @(k1,k2) (((M(k1,k2,2,2)-M(k1,k2,1,1))./2)-SQ(k1,k2))./sqrt(M(k1,k2,1,2).*M(k1,k2,2,1)+(((M(k1,k2,2,2)-M(k1,k2,1,1))./2)-SQ(k1,k2)).^2);
xp= @(k1,k2) M(k1,k2,1,2)./sqrt(M(k1,k2,1,2).*M(k1,k2,2,1)+(((M(k1,k2,2,2)-M(k1,k2,1,1))./2)+SQ(k1,k2)).^2);
xm= @(k1,k2) M(k1,k2,1,2)./sqrt(M(k1,k2,1,2).*M(k1,k2,2,1)+(((M(k1,k2,2,2)-M(k1,k2,1,1))./2)-SQ(k1,k2)).^2);

Lm= @(k1,k2) (xm(k1,k2).*transpose(k1)+ym(k1,k2).*transpose(k2))...
              .*conj((xm(k1,k2).*transpose(k1)+ym(k1,k2).*transpose(k2)))...
              ./transpose(k1.^2+k2.^2);
Lp= @(k1,k2) 1-Lm(k1,k2);
end



