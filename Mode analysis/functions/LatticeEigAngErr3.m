function [S,Err,Errext]=LatticeEigAngErr3(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x,ang,l0,frac,Weight)
Diff=frac.*x;

S0=LatticeFitAng3(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x,ang,l0,0);

d00=sqrt(LatticeFitAng3(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x,ang,l0,Weight));
d10=sqrt(LatticeFitAng3(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x+[Diff(1),0],ang,l0,Weight));
dm0=sqrt(LatticeFitAng3(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x-[Diff(1),0],ang,l0,Weight));
d01=sqrt(LatticeFitAng3(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x+[0,Diff(2)],ang,l0,Weight));
d0m=sqrt(LatticeFitAng3(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x-[0,Diff(2)],ang,l0,Weight));
d11=sqrt(LatticeFitAng3(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x+[Diff(1),Diff(2)],ang,l0,Weight));
dmm=sqrt(LatticeFitAng3(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x-[Diff(1),Diff(2)],ang,l0,Weight));
d1m=sqrt(LatticeFitAng3(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x+[Diff(1),-Diff(2)],ang,l0,Weight));
dm1=sqrt(LatticeFitAng3(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x+[-Diff(1),Diff(2)],ang,l0,Weight));

S=[(d10-dm0)/(2*Diff(1)),(d01-d0m)/(2*Diff(2))];

Sxx=(d10+dm0-2*d00)/Diff(1)^2;
Sxy=(d11+dmm-d1m-dm1)/(4*Diff(1)*Diff(2));
Syy=(d01+d0m-2*d00)/Diff(2)^2;

[lp,lm]=LatticeEigAng(x,ang,l0);
[lp10,lm10]=LatticeEigAng(x+[Diff(1),0],ang,l0);
[lpm0,lmm0]=LatticeEigAng(x-[Diff(1),0],ang,l0);
[lp01,lm01]=LatticeEigAng(x+[0,Diff(2)],ang,l0);
[lp0m,lm0m]=LatticeEigAng(x-[0,Diff(2)],ang,l0);

dfp1= @(k1,k2) (sqrt(lp10(k1,k2))-sqrt(lpm0(k1,k2)))/(2*Diff(1));
dfp2= @(k1,k2) (sqrt(lp01(k1,k2))-sqrt(lp0m(k1,k2)))/(2*Diff(2));
dfm1= @(k1,k2) (sqrt(lm10(k1,k2))-sqrt(lmm0(k1,k2)))/(2*Diff(1));
dfm2= @(k1,k2) (sqrt(lm01(k1,k2))-sqrt(lm0m(k1,k2)))/(2*Diff(2));


Err=[0.0 0.0];
Errext=[0.0 0.0];
Tot=sum(1-Echeck);
for i=1:length(K)
    if Echeck(i)==0
        %D1=real(((sqrt(lp(K(i,1),K(i,2)))-Wp(i))/real(Wpperr(i)))^2+((sqrt(lm(K(i,1),K(i,2)))-Wm(i))/real(Wmmerr(i)))^2);
        %D2=real(((sqrt(lm(K(i,1),K(i,2)))-Wp(i))/real(Wpperr(i)))^2+((sqrt(lp(K(i,1),K(i,2)))-Wm(i))/real(Wmmerr(i)))^2);
        if Weight==1
            dm1dwp=(d00*2*Tot*(Sxx*Syy-Sxy^2)*(Wpperr(i)*Wmmerr(i)-Wpmerr(i).^2)).^-1*(Wpperr(i)*(Syy*dfp1(K(i,1),K(i,2))-Sxy*dfp2(K(i,1),K(i,2)))+Wpmerr(i)*(Syy*dfm1(K(i,1),K(i,2))-Sxy*dfm2(K(i,1),K(i,2))));
            dm1dwm=(d00*2*Tot*(Sxx*Syy-Sxy^2)*(Wpperr(i)*Wmmerr(i)-Wpmerr(i).^2)).^-1*(Wpmerr(i)*(Syy*dfp1(K(i,1),K(i,2))-Sxy*dfp2(K(i,1),K(i,2)))+Wmmerr(i)*(Syy*dfm1(K(i,1),K(i,2))-Sxy*dfm2(K(i,1),K(i,2))));
            dm2dwp=(d00*2*Tot*(Sxx*Syy-Sxy^2)*(Wpperr(i)*Wmmerr(i)-Wpmerr(i).^2)).^-1*(Wpperr(i)*(Sxx*dfp2(K(i,1),K(i,2))-Sxy*dfp1(K(i,1),K(i,2)))+Wpmerr(i)*(Sxx*dfm2(K(i,1),K(i,2))-Sxy*dfm1(K(i,1),K(i,2))));
            dm2dwm=(d00*2*Tot*(Sxx*Syy-Sxy^2)*(Wpperr(i)*Wmmerr(i)-Wpmerr(i).^2)).^-1*(Wpmerr(i)*(Sxx*dfp2(K(i,1),K(i,2))-Sxy*dfp1(K(i,1),K(i,2)))+Wmmerr(i)*(Sxx*dfm2(K(i,1),K(i,2))-Sxy*dfm1(K(i,1),K(i,2))));
            %recall that plus and minus switch in the data vs calculated
            %since there is an inverse.
            Err(1)=Err(1)+dm1dwp.*Wmmerr(i).*dm1dwp+2.*dm1dwp.*Wpmerr(i).*dm1dwm+dm1dwm*Wpperr(i)*dm1dwm;
            Err(2)=Err(2)+dm2dwp.*Wmmerr(i).*dm2dwp+2.*dm2dwp.*Wpmerr(i).*dm2dwm+dm2dwm*Wpperr(i)*dm2dwm;
            %Errext(1)=Errext(1)+(dm1dwp.*(Wm(i)-sqrt(lp(K(i,1),K(i,2)))).^2.*dm1dwp+dm1dwm*(Wp(i)-sqrt(lm(K(i,1),K(i,2)))).^2*dm1dwm)./(2*N);
            %Errext(2)=Errext(2)+(dm2dwp.*(Wm(i)-sqrt(lp(K(i,1),K(i,2)))).^2.*dm2dwp+dm2dwm*(Wp(i)-sqrt(lm(K(i,1),K(i,2)))).^2*dm2dwm)./(2*N);
            Errext(1)=Errext(1)+(dm1dwp.*dm1dwp+dm1dwm.*dm1dwm)*S0;
            Errext(2)=Errext(2)+(dm2dwp.*dm2dwp+dm2dwm.*dm2dwm)*S0;
            %Err(1)=Err(1)+(d00*Tot*2*(Sxx*Syy-Sxy^2)/2).^-2*((Syy*dfp1(K(i,1),K(i,2))-Sxy*dfp2(K(i,1),K(i,2))).*Wmmerr(i).*(Syy*dfp1(K(i,1),K(i,2))-Sxy*dfp2(K(i,1),K(i,2)))...
            %    -2*(Syy*dfp1(K(i,1),K(i,2))-Sxy*dfp2(K(i,1),K(i,2))).*Wpmerr(i).*(Syy*dfm1(K(i,1),K(i,2))-Sxy*dfm2(K(i,1),K(i,2)))...
            %    +(Syy*dfm1(K(i,1),K(i,2))-Sxy*dfm2(K(i,1),K(i,2))).*Wpperr(i).*(Syy*dfm1(K(i,1),K(i,2))-Sxy*dfm2(K(i,1),K(i,2))))/(Wpperr(i)*Wmmerr(i)-Wpmerr(i)^2);
            %Err(2)=Err(2)+(d00*Tot*2*(Sxx*Syy-Sxy^2)/2).^-2*((Sxx*dfp2(K(i,1),K(i,2))-Sxy*dfp1(K(i,1),K(i,2))).*Wmmerr(i).*(Sxx*dfp2(K(i,1),K(i,2))-Sxy*dfp1(K(i,1),K(i,2)))...
            %    -2*(Sxx*dfp2(K(i,1),K(i,2))-Sxy*dfp1(K(i,1),K(i,2))).*Wpmerr(i).*(Sxx*dfm2(K(i,1),K(i,2))-Sxy*dfm1(K(i,1),K(i,2)))...
            %    +(Sxx*dfm2(K(i,1),K(i,2))-Sxy*dfm1(K(i,1),K(i,2))).*Wpperr(i).*(Sxx*dfm2(K(i,1),K(i,2))-Sxy*dfm1(K(i,1),K(i,2))))/(Wpperr(i)*Wmmerr(i)-Wpmerr(i)^2);
        else
            dm1dwp=(d00*2*Tot*(Sxx*Syy-Sxy^2)).^-1*(Syy*dfp1(K(i,1),K(i,2))-Sxy*dfp2(K(i,1),K(i,2)));
            dm1dwm=(d00*2*Tot*(Sxx*Syy-Sxy^2)).^-1*(Syy*dfm1(K(i,1),K(i,2))-Sxy*dfm2(K(i,1),K(i,2)));
            dm2dwp=(d00*2*Tot*(Sxx*Syy-Sxy^2)).^-1*(Sxx*dfp2(K(i,1),K(i,2))-Sxy*dfp1(K(i,1),K(i,2)));
            dm2dwm=(d00*2*Tot*(Sxx*Syy-Sxy^2)).^-1*(Sxx*dfm2(K(i,1),K(i,2))-Sxy*dfm1(K(i,1),K(i,2)));
            Err(1)=Err(1)+dm1dwp.*Wmmerr(i).*dm1dwp+2.*dm1dwp.*Wpmerr(i).*dm1dwm+dm1dwm*Wpperr(i)*dm1dwm;
            Err(2)=Err(2)+dm2dwp.*Wmmerr(i).*dm2dwp+2.*dm2dwp.*Wpmerr(i).*dm2dwm+dm2dwm*Wpperr(i)*dm2dwm;
            Errext(1)=Errext(1)+(dm1dwp.*dm1dwp+dm1dwm.*dm1dwm)*S0;
            Errext(2)=Errext(2)+(dm2dwp.*dm2dwp+dm2dwm.*dm2dwm)*S0;
            %Errext(1)=Errext(1)+(dm1dwp.*(Wm(i)-sqrt(lp(K(i,1),K(i,2)))).^2.*dm1dwp+dm1dwm*(Wp(i)-sqrt(lm(K(i,1),K(i,2)))).^2*dm1dwm);
            %Errext(2)=Errext(2)+(dm2dwp.*(Wm(i)-sqrt(lp(K(i,1),K(i,2)))).^2.*dm2dwp+dm2dwm*(Wp(i)-sqrt(lm(K(i,1),K(i,2)))).^2*dm2dwm);
        end
    end
end
Err=sqrt(Err);%d00*sqrt(Err);




