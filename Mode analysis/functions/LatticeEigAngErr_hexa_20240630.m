function [S,Err,Errext]=LatticeEigAngErr_hexa_20240630(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x,l0,frac,Weight)

[lm,lp]=LatticeEig_hexa(x,l0);

% Sum of normalized, squared residuals, chi^2(k)
chi2 = ((Wp - lp(K(:,1),K(:,2))')./Wpperr).^2 + ((Wm - lm(K(:,1),K(:,2))')./Wmmerr).^2;
chi2(Echeck'==1)=0;
chi2 = sum(chi2);

% Number of data points
N = length(Wp);

% 1st derivative of residuals, with analytical solution.
[lm1,lp1]=LatticeEig_hexa(1,l0);
C = (lp1(K(:,1),K(:,2))'./Wpperr).^2 + (lm1(K(:,1),K(:,2))'./Wmmerr).^2;
C(Echeck'==1) = 0;
C = sum(C,"all");

% Full error
Err = chi2 / (2*N-1) / C;
Err = sqrt(Err);%d00*sqrt(Err);

S = 0;
Errext = 0;

end


% Diff=frac.*x;
% 
% S0=LatticeFit_hexa(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x,l0,0);
% 
% chi2=LatticeFit_hexa(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x,l0,Weight);
% d10=LatticeFit_hexa(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x+Diff,l0,Weight);
% dm0=LatticeFit_hexa(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x-Diff,l0,Weight);
% 
% S=(d10-dm0)/(2*Diff); % 1st derivative
% 
% Sxx=(d10+dm0-2*chi2)/Diff^2;
% 
% 
% [lp10,lm10]=LatticeEig_hexa(x+Diff,l0);
% [lpm0,lmm0]=LatticeEig_hexa(x-Diff,l0);
% dfp1= @(k1,k2) (sqrt(lp10(k1,k2))-sqrt(lpm0(k1,k2)))/(2*Diff(1));
% dfm1= @(k1,k2) (sqrt(lm10(k1,k2))-sqrt(lmm0(k1,k2)))/(2*Diff(1));
% C = dfp1(K(:,1), K(:,2))' ./ Wpperr + dfm1(K(:,1), K(:,2))' ./ Wmmerr;
% C(Echeck'==1) = 0;
% C = sum(C.^2,"all");
% 
% 
% Err = chi2.^2 / (2*N-1) / C;
% Err = sqrt(Err);%d00*sqrt(Err);
% 
% Errext = 0;


% 
% 
% Err=[0.0];
% Errext=[0.0];
% Tot=sum(1-Echeck);
% for i=1:length(K)
%     if Echeck(i)==0
%         if Weight==1
%             dm1dwp=(d00*2*Tot*(Sxx*Syy-Sxy^2)*(Wpperr(i)*Wmmerr(i)-Wpmerr(i).^2)).^-1*(Wpperr(i)*(Syy*dfp1(K(i,1),K(i,2))-Sxy*dfp2(K(i,1),K(i,2)))+Wpmerr(i)*(Syy*dfm1(K(i,1),K(i,2))-Sxy*dfm2(K(i,1),K(i,2))));
%             dm1dwm=(d00*2*Tot*(Sxx*Syy-Sxy^2)*(Wpperr(i)*Wmmerr(i)-Wpmerr(i).^2)).^-1*(Wpmerr(i)*(Syy*dfp1(K(i,1),K(i,2))-Sxy*dfp2(K(i,1),K(i,2)))+Wmmerr(i)*(Syy*dfm1(K(i,1),K(i,2))-Sxy*dfm2(K(i,1),K(i,2))));
%             dm2dwp=(d00*2*Tot*(Sxx*Syy-Sxy^2)*(Wpperr(i)*Wmmerr(i)-Wpmerr(i).^2)).^-1*(Wpperr(i)*(Sxx*dfp2(K(i,1),K(i,2))-Sxy*dfp1(K(i,1),K(i,2)))+Wpmerr(i)*(Sxx*dfm2(K(i,1),K(i,2))-Sxy*dfm1(K(i,1),K(i,2))));
%             dm2dwm=(d00*2*Tot*(Sxx*Syy-Sxy^2)*(Wpperr(i)*Wmmerr(i)-Wpmerr(i).^2)).^-1*(Wpmerr(i)*(Sxx*dfp2(K(i,1),K(i,2))-Sxy*dfp1(K(i,1),K(i,2)))+Wmmerr(i)*(Sxx*dfm2(K(i,1),K(i,2))-Sxy*dfm1(K(i,1),K(i,2))));
%             %recall that plus and minus switch in the data vs calculated
%             %since there is an inverse.
%             Err(1)=Err(1)+dm1dwp.*Wmmerr(i).*dm1dwp+2.*dm1dwp.*Wpmerr(i).*dm1dwm+dm1dwm*Wpperr(i)*dm1dwm;
%             Err(2)=Err(2)+dm2dwp.*Wmmerr(i).*dm2dwp+2.*dm2dwp.*Wpmerr(i).*dm2dwm+dm2dwm*Wpperr(i)*dm2dwm;
% 
%             Errext(1)=Errext(1)+(dm1dwp.*dm1dwp+dm1dwm.*dm1dwm)*S0;
%             Errext(2)=Errext(2)+(dm2dwp.*dm2dwp+dm2dwm.*dm2dwm)*S0;
%         else
%             % dm1dwp=(d00*2*Tot*(Sxx^2)).^-1;
%             % dm1dwm=(d00*2*Tot*(Sxx^2)).^-1;
%             dm1dwp=(d00*2*Tot*(Sxx^2)).^-1*(Sxx*dfp1(K(i,1),K(i,2)));
%             dm1dwm=(d00*2*Tot*(Sxx^2)).^-1*(Sxx*dfm1(K(i,1),K(i,2)));
%             Err(1)=Err(1)+dm1dwp.*Wmmerr(i).*dm1dwp+2.*dm1dwp.*Wpmerr(i).*dm1dwm+dm1dwm*Wpperr(i)*dm1dwm;
%             Errext(1)=Errext(1)+(dm1dwp.*dm1dwp+dm1dwm.*dm1dwm)*S0;
%         end
%     end
% end
% Err=sqrt(Err);%d00*sqrt(Err);




