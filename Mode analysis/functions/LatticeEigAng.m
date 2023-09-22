function [lp,lm]=LatticeEigAng(x,ang,l0)
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


%Matrix
%a=ang
% M=@(kx,ky,a,b) 2*k*(a1(a)*a1(b)+a2(a)*a2(b))-2.*k.*a1(a).*a1(b).*cos(kx*A1(1)+ky*A1(2))-2.*k.*a2(a).*a2(b).*cos(kx*A2(1)+ky*A2(2))+Q*(a1p(a)-a2p(a))*(a1p(b)-a2p(b))+Q*(a1p(a))*(a1p(b))+Q*(a2p(a))*(a2p(b))+Q*(a1p(a)-a2p(a))*(-a1p(b))*exp(j*(kx*A1(1)+ky*A1(2)))+Q*(a1p(b)-a2p(b))*(-a1p(a))*exp(-j*(kx*A1(1)+ky*A1(2)))+Q*(a1p(a)-a2p(a))*(-a2p(b))*exp(j*(kx*A2(1)+ky*A2(2)))+Q*(a1p(b)-a2p(b))*(-a2p(a))*exp(-j*(kx*A2(1)+ky*A2(2)))+Q*(-a1p(a))*(a2p(b))*exp(-j*(kx*B2(1)+ky*B2(2)))+Q*(-a1p(b))*(a2p(a))*exp(j*(kx*B2(1)+ky*B2(2)))+Q*(a2p(a)+a1p(a))*(a2p(b)+a1p(b))+Q*(a1p(a))*(a1p(b))+Q*(a2p(a))*(a2p(b))+Q*(a2p(a)+a1p(a))*(-a2p(b))*exp(j*(kx*A2(1)+ky*A2(2)))+Q*(a2p(b)+a1p(b))*(-a2p(a))*exp(-j*(kx*A2(1)+ky*A2(2)))+Q*(a2p(a)+a1p(a))*(a1p(b))*exp(-j*(kx*A1(1)+ky*A1(2)))+Q*(a2p(b)+a1p(b))*(a1p(a))*exp(j*(kx*A1(1)+ky*A1(2)))+Q*(-a2p(a))*(-a1p(b))*exp(-j*(kx*B1(1)+ky*B1(2)))+Q*(-a2p(b))*(-a1p(a))*exp(j*(kx*B1(1)+ky*B1(2)))+Q*(a1p(a)-a2p(a))*(a1p(b)-a2p(b))+Q*(a1p(a))*(a1p(b))+Q*(a2p(a))*(a2p(b))+Q*(a1p(a)-a2p(a))*(-a1p(b))*exp(-j*(kx*A1(1)+ky*A1(2)))+Q*(a1p(b)-a2p(b))*(-a1p(a))*exp(j*(kx*A1(1)+ky*A1(2)))+Q*(a1p(a)-a2p(a))*(-a2p(b))*exp(-j*(kx*A2(1)+ky*A2(2)))+Q*(a1p(b)-a2p(b))*(-a2p(a))*exp(j*(kx*A2(1)+ky*A2(2)))+Q*(-a1p(a))*(a2p(b))*exp(j*(kx*B2(1)+ky*B2(2)))+Q*(-a1p(b))*(a2p(a))*exp(-j*(kx*B2(1)+ky*B2(2)))+Q*(a2p(a)+a1p(a))*(a2p(b)+a1p(b))+Q*(a1p(a))*(a1p(b))+Q*(a2p(a))*(a2p(b))+Q*(a2p(a)+a1p(a))*(-a2p(b))*exp(-j*(kx*A2(1)+ky*A2(2)))+Q*(a2p(b)+a1p(b))*(-a2p(a))*exp(j*(kx*A2(1)+ky*A2(2)))+Q*(a2p(a)+a1p(a))*(a1p(b))*exp(j*(kx*A1(1)+ky*A1(2)))+Q*(a2p(b)+a1p(b))*(a1p(a))*exp(-j*(kx*A1(1)+ky*A1(2)))+Q*(-a2p(a))*(-a1p(b))*exp(j*(kx*B1(1)+ky*B1(2)))+Q*(-a2p(b))*(-a1p(a))*exp(-j*(kx*B1(1)+ky*B1(2)));
% M=@(kx,ky,a,b) 2*k*(a1(a)*a1(b)+a2(a)*a2(b))-2.*k.*a1(a).*a1(b).*cos(kx*A1(1)+ky*A1(2))-2.*k.*a2(a).*a2(b).*cos(kx*A2(1)+ky*A2(2))+Q*(a1p(a)-a2p(a))*(a1p(b)-a2p(b))+Q*(a1p(a))*(a1p(b))+Q*(a2p(a))*(a2p(b))+Q*(a1p(a)-a2p(a))*(-a1p(b))*exp(j*(kx*A1(1)+ky*A1(2)))+Q*(a1p(b)-a2p(b))*(-a1p(a))*exp(-j*(kx*A1(1)+ky*A1(2)))+Q*(a1p(a)-a2p(a))*(a2p(b))*exp(j*(kx*A2(1)+ky*A2(2)))+Q*(a1p(b)-a2p(b))*(a2p(a))*exp(-j*(kx*A2(1)+ky*A2(2)))+Q*(-a1p(a))*(a2p(b))*exp(-j*(kx*B2(1)+ky*B2(2)))+Q*(-a1p(b))*(a2p(a))*exp(j*(kx*B2(1)+ky*B2(2)))+Q*(a2p(a)+a1p(a))*(a2p(b)+a1p(b))+Q*(a1p(a))*(a1p(b))+Q*(a2p(a))*(a2p(b))+Q*(a2p(a)+a1p(a))*(-a2p(b))*exp(j*(kx*A2(1)+ky*A2(2)))+Q*(a2p(b)+a1p(b))*(-a2p(a))*exp(-j*(kx*A2(1)+ky*A2(2)))+Q*(a2p(a)+a1p(a))*(-a1p(b))*exp(-j*(kx*A1(1)+ky*A1(2)))+Q*(a2p(b)+a1p(b))*(-a1p(a))*exp(j*(kx*A1(1)+ky*A1(2)))+Q*(-a2p(a))*(-a1p(b))*exp(-j*(kx*B1(1)+ky*B1(2)))+Q*(-a2p(b))*(-a1p(a))*exp(j*(kx*B1(1)+ky*B1(2)))+Q*(a1p(a)-a2p(a))*(a1p(b)-a2p(b))+Q*(a1p(a))*(a1p(b))+Q*(a2p(a))*(a2p(b))+Q*(a1p(a)-a2p(a))*(-a1p(b))*exp(-j*(kx*A1(1)+ky*A1(2)))+Q*(a1p(b)-a2p(b))*(-a1p(a))*exp(j*(kx*A1(1)+ky*A1(2)))+Q*(a1p(a)-a2p(a))*(a2p(b))*exp(-j*(kx*A2(1)+ky*A2(2)))+Q*(a1p(b)-a2p(b))*(a2p(a))*exp(j*(kx*A2(1)+ky*A2(2)))+Q*(-a1p(a))*(a2p(b))*exp(j*(kx*B2(1)+ky*B2(2)))+Q*(-a1p(b))*(a2p(a))*exp(-j*(kx*B2(1)+ky*B2(2)))+Q*(a2p(a)+a1p(a))*(a2p(b)+a1p(b))+Q*(a1p(a))*(a1p(b))+Q*(a2p(a))*(a2p(b))+Q*(a2p(a)+a1p(a))*(-a2p(b))*exp(-j*(kx*A2(1)+ky*A2(2)))+Q*(a2p(b)+a1p(b))*(-a2p(a))*exp(j*(kx*A2(1)+ky*A2(2)))+Q*(a2p(a)+a1p(a))*(-a1p(b))*exp(j*(kx*A1(1)+ky*A1(2)))+Q*(a2p(b)+a1p(b))*(-a1p(a))*exp(-j*(kx*A1(1)+ky*A1(2)))+Q*(-a2p(a))*(-a1p(b))*exp(j*(kx*B1(1)+ky*B1(2)))+Q*(-a2p(b))*(-a1p(a))*exp(-j*(kx*B1(1)+ky*B1(2)));
% 
% A=@(kx,ky) M(kx,ky,1,1);
% B1=@(kx,ky) M(kx,ky,1,2);
% B2=@(kx,ky) M(kx,ky,2,1);
% C=@(kx,ky) M(kx,ky,2,2);
% %Eigenvalues
% lp= @(k1,k2) 0.5.*((A(k1,k2)+C(k1,k2))+sqrt((A(k1,k2)+C(k1,k2)).^2-4.*(A(k1,k2).*C(k1,k2) - B1(k1,k2).*B2(k1,k2) ) ));
% lm= @(k1,k2) 0.5.*((A(k1,k2)+C(k1,k2))-sqrt((A(k1,k2)+C(k1,k2)).^2-4.*(A(k1,k2).*C(k1,k2) - B1(k1,k2).*B2(k1,k2) ) ));

M=@(kx,ky,a,b) 2*k*(a1(a)*a1(b).*(1-cos(kx*A1(1)+ky*A1(2)))+a2(a)*a2(b).*(1-cos(kx*A2(1)+ky*A2(2)))) ...
    +16*x(2)*(2*ang-pi)^2/(l0^2*sin(ang)^2).*(cos(kx*A1(1)+ky*A1(2))-1).*(cos(ang)*(a1(a)*a2(b)+a2(a)*a1(b))-cos(ang)^2*a1(a)*a1(b)-a2(a)*a2(b))...
    +16*x(2)*(2*ang-pi)^2/(l0^2*sin(ang)^2).*(cos(kx*A2(1)+ky*A2(2))-1).*(cos(ang)*(a1(a)*a2(b)+a2(a)*a1(b))-cos(ang)^2*a2(a)*a2(b)-a1(a)*a1(b))...
    +4*x(2)*(2*ang-pi)^2/(l0^2*sin(ang)^2).*(cos(kx*(A1(1)+A2(1))+ky*(A1(2)+A2(2)))-cos(kx*(A1(1)-A2(1))+ky*(A1(2)-A2(2)))).*(2*cos(ang)*(a2(a)*a2(b)+a1(a)*a1(b))-(1+cos(ang)^2)*(a1(a)*a2(b)+a2(a)*a1(b)));

lp= @(k1,k2) 0.5*((M(k1,k2,1,1)+M(k1,k2,2,2))+sqrt((M(k1,k2,1,1)+M(k1,k2,2,2)).^2-4*(M(k1,k2,1,1).*M(k1,k2,2,2) - M(k1,k2,2,1).*M(k1,k2,1,2) ) ));
lm= @(k1,k2) 0.5*((M(k1,k2,1,1)+M(k1,k2,2,2))-sqrt((M(k1,k2,1,1)+M(k1,k2,2,2)).^2-4*(M(k1,k2,1,1).*M(k1,k2,2,2) - M(k1,k2,2,1).*M(k1,k2,1,2) ) ));






