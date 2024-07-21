function [d]=LatticeFitAng3check(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x,a,l0,Weight)
[lp,lm]=LatticeEigAng(x,a,l0);
%fitting
d=0.0;
dm=0.0;
dp=0.0;
de=0.0;
Tot=sum(1-Echeck);
for i=1:length(K)
    if Echeck(i)==0
        if Weight==1%standard weighted least square (Error not computed as of yet)
            D=((sqrt(lm(K(i,1),K(i,2)))-Wp(i))*Wmmerr(i)*(sqrt(lm(K(i,1),K(i,2)))-Wp(i))...
            -2*(sqrt(lm(K(i,1),K(i,2)))-Wp(i))*Wpmerr(i)*(sqrt(lp(K(i,1),K(i,2)))-Wm(i))...
            +(sqrt(lp(K(i,1),K(i,2)))-Wm(i))*Wpperr(i)*(sqrt(lp(K(i,1),K(i,2)))-Wm(i)))/(Wpperr(i)*Wmmerr(i)-Wpmerr(i)^2);
        else %unweighted
            Dm=(sqrt(lm(K(i,1),K(i,2)))-Wp(i)).^2;
            Dp=(sqrt(lp(K(i,1),K(i,2)))-Wm(i)).^2;
            De=(sqrt(lm(K(i,1),K(i,2)))-Wp(i))*Wpperr(i)*(sqrt(lm(K(i,1),K(i,2)))-Wp(i))...
            +2*(sqrt(lm(K(i,1),K(i,2)))-Wp(i))*Wpmerr(i)*(sqrt(lp(K(i,1),K(i,2)))-Wm(i))...
            +(sqrt(lp(K(i,1),K(i,2)))-Wm(i))*Wmmerr(i)*(sqrt(lp(K(i,1),K(i,2)))-Wm(i));
        end

        if Weight==1
            d = d+D/Tot/2;
        else
        d=d+(Dm+Dp)/Tot/2;
        end
        % dm=dm+Dm/Tot/2;
        % dp=dp+Dp/Tot/2;
        % de=de+De/(Tot*2)^2;
        % if Dm+Dp<0
        %     disp(x)
        %     disp([lm(K(i,1),K(i,2)),lp(K(i,1),K(i,2)),D/Tot/2])
        % end
    end
end
d=sqrt(d);
% dm=sqrt(dm);
% dp=sqrt(dp);
% de=sqrt(de)/d;

