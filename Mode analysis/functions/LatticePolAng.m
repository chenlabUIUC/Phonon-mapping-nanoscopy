% function [Lp,Lm,xp,yp,xm,ym]=LatticePolAng(x,ang,l0)
function [Lp,Lm]=LatticePolAng(x,ang,l0)
%x=[k,a,l0,E]
%Secondary Constants
x=abs(x);
k=x(1);
%Angular:
Q=2*x(2)*(2*ang-pi)^2/l0^2;
%j=sqrt(-1);
%Vectors
A1=[l0,0];
A2=[l0*cos(ang),l0*sin(ang)];
a1=[1,0];
a2=[cos(ang),sin(ang)];
%diagonals:
B1=A1+A2;
B2=A1-A2;
% b1=B1./(2*l0*cos(ang/2));
% b2=B2./(2*l0*sin(ang/2));
%perpendicular vectors
a1p=[0,1];
a2p=[-sin(ang),cos(ang)];

M=@(kx,ky,a,b) 2*k*(a1(a)*a1(b).*(1-cos(kx*A1(1)+ky*A1(2)))+a2(a)*a2(b).*(1-cos(kx*A2(1)+ky*A2(2)))) ...
    +16*x(2)*(2*ang-pi)^2/(l0^2*sin(ang)^2).*(cos(kx*A1(1)+ky*A1(2))-1).*(cos(ang)*(a1(a)*a2(b)+a2(a)*a1(b))-cos(ang)^2*a1(a)*a1(b)-a2(a)*a2(b))...
    +16*x(2)*(2*ang-pi)^2/(l0^2*sin(ang)^2).*(cos(kx*A2(1)+ky*A2(2))-1).*(cos(ang)*(a1(a)*a2(b)+a2(a)*a1(b))-cos(ang)^2*a2(a)*a2(b)-a1(a)*a1(b))...
    +4*x(2)*(2*ang-pi)^2/(l0^2*sin(ang)^2).*(cos(kx*(A1(1)+A2(1))+ky*(A1(2)+A2(2)))-cos(kx*(A1(1)-A2(1))+ky*(A1(2)-A2(2)))).*(2*cos(ang)*(a2(a)*a2(b)+a1(a)*a1(b))-(1+cos(ang)^2)*(a1(a)*a2(b)+a2(a)*a1(b)));

%lp= @(k1,k2) 0.5*((M(k1,k2,1,1)+M(k1,k2,2,2))+sqrt((M(k1,k2,1,1)+M(k1,k2,2,2)).^2-4*(M(k1,k2,1,1).*M(k1,k2,2,2) - M(k1,k2,2,1).*M(k1,k2,1,2) ) ));
%lm= @(k1,k2) 0.5*((M(k1,k2,1,1)+M(k1,k2,2,2))-sqrt((M(k1,k2,1,1)+M(k1,k2,2,2)).^2-4*(M(k1,k2,1,1).*M(k1,k2,2,2) - M(k1,k2,2,1).*M(k1,k2,1,2) ) ));

SQ= @(k1,k2) sqrt(((M(k1,k2,1,1)+M(k1,k2,2,2))./2).^2-(M(k1,k2,1,1).*M(k1,k2,2,2)-M(k1,k2,1,2).*M(k1,k2,2,1)));

yp= @(k1,k2) (((M(k1,k2,2,2)-M(k1,k2,1,1))./2)+SQ(k1,k2))./sqrt(M(k1,k2,1,2).*M(k1,k2,2,1)+(((M(k1,k2,2,2)-M(k1,k2,1,1))./2)+SQ(k1,k2)).^2);
ym= @(k1,k2) (((M(k1,k2,2,2)-M(k1,k2,1,1))./2)-SQ(k1,k2))./sqrt(M(k1,k2,1,2).*M(k1,k2,2,1)+(((M(k1,k2,2,2)-M(k1,k2,1,1))./2)-SQ(k1,k2)).^2);
xp= @(k1,k2) M(k1,k2,1,2)./sqrt(M(k1,k2,1,2).*M(k1,k2,2,1)+(((M(k1,k2,2,2)-M(k1,k2,1,1))./2)+SQ(k1,k2)).^2);
xm= @(k1,k2) M(k1,k2,1,2)./sqrt(M(k1,k2,1,2).*M(k1,k2,2,1)+(((M(k1,k2,2,2)-M(k1,k2,1,1))./2)-SQ(k1,k2)).^2);

Lm= @(k1,k2) (xm(k1,k2).*transpose(k1)+ym(k1,k2).*transpose(k2))...
              .*conj((xm(k1,k2).*transpose(k1)+ym(k1,k2).*transpose(k2)))...
              ./transpose(k1.^2+k2.^2);
Lp= @(k1,k2) 1-Lm(k1,k2);
% Lm= @(k1,k2) (xm(k1,k2).*conj(xm(k1,k2)).*transpose(k1).^2+ym(k1,k2).^2.*transpose(k2).^2)./transpose(k1.^2+k2.^2);
% Lp= @(k1,k2) 1-Lm(k1,k2)%(xp(k1,k2).*conj(xp(k1,k2)).*transpose(k1).^2+yp(k1,k2).^2.*transpose(k2).^2)./transpose(k1.^2+k2.^2);




