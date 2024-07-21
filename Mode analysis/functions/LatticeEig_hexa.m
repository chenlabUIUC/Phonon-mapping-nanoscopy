function [lp,lm]=LatticeEig_hexa(x,l0)
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
% M=@(kx,ky,a,b) 2*k*(a1(a)*a1(b).*(1-cos(kx*A1(1)+ky*A1(2)))+a2(a)*a2(b).*(1-cos(kx*A2(1)+ky*A2(2)))) ...
%     +16*x(2)*(2*ang-pi)^2/(l0^2*sin(ang)^2).*(cos(kx*A1(1)+ky*A1(2))-1).*(cos(ang)*(a1(a)*a2(b)+a2(a)*a1(b))-cos(ang)^2*a1(a)*a1(b)-a2(a)*a2(b))...
%     +16*x(2)*(2*ang-pi)^2/(l0^2*sin(ang)^2).*(cos(kx*A2(1)+ky*A2(2))-1).*(cos(ang)*(a1(a)*a2(b)+a2(a)*a1(b))-cos(ang)^2*a2(a)*a2(b)-a1(a)*a1(b))...
%     +4*x(2)*(2*ang-pi)^2/(l0^2*sin(ang)^2).*(cos(kx*(A1(1)+A2(1))+ky*(A1(2)+A2(2)))-cos(kx*(A1(1)-A2(1))+ky*(A1(2)-A2(2)))).*(2*cos(ang)*(a2(a)*a2(b)+a1(a)*a1(b))-(1+cos(ang)^2)*(a1(a)*a2(b)+a2(a)*a1(b)));

M=@(kx,ky,a,b) 2*k*(a1(a)*a1(b).*(1-cos(kx*A1(1)+ky*A1(2)))+...
                    a2(a)*a2(b).*(1-cos(kx*A2(1)+ky*A2(2)))+...
                    a3(a)*a3(b).*(1-cos(kx*A3(1)+ky*A3(2))));

lp= @(k1,k2) 0.5*((M(k1,k2,1,1)+M(k1,k2,2,2))+sqrt((M(k1,k2,1,1)+M(k1,k2,2,2)).^2-4*(M(k1,k2,1,1).*M(k1,k2,2,2) - M(k1,k2,2,1).*M(k1,k2,1,2) ) ));
lm= @(k1,k2) 0.5*((M(k1,k2,1,1)+M(k1,k2,2,2))-sqrt((M(k1,k2,1,1)+M(k1,k2,2,2)).^2-4*(M(k1,k2,1,1).*M(k1,k2,2,2) - M(k1,k2,2,1).*M(k1,k2,1,2) ) ));






