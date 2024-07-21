function d=LatticeFitAng4(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x,a,l0,Weight,only_lower)
[lp,lm]=LatticeEigAng(x,a,l0);
%
if only_lower==1
    Wm = zeros(size(Wm));
    for i = 1:size(K,1)
        Wm(i) = sqrt(lp(K(i,1),K(i,2)));
    end
end

%fitting
d=0.0;
Tot=sum(1-Echeck);
for i=1:length(K)
    if Echeck(i)==0
        if Weight==1%standard weighted least square (more correlation than is currently accounted for?)
            D=((sqrt(lm(K(i,1),K(i,2)))-Wp(i))*Wmmerr(i)*(sqrt(lm(K(i,1),K(i,2)))-Wp(i))...
            -2*(sqrt(lm(K(i,1),K(i,2)))-Wp(i))*Wpmerr(i)*(sqrt(lp(K(i,1),K(i,2)))-Wm(i))...
            +(sqrt(lp(K(i,1),K(i,2)))-Wm(i))*Wpperr(i)*(sqrt(lp(K(i,1),K(i,2)))-Wm(i)))/(Wpperr(i)*Wmmerr(i)-Wpmerr(i)^2);
        else %unweighted
            D=((sqrt(lm(K(i,1),K(i,2)))-Wp(i)).^2+(sqrt(lp(K(i,1),K(i,2)))-Wm(i)).^2);
        end
            d=d+D/Tot/2;
        if D<0
            disp(x)
            disp([lm(K(i,1),K(i,2)),lp(K(i,1),K(i,2)),D/Tot/2])
        end
    end
end

