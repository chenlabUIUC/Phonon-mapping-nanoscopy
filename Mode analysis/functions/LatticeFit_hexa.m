function d = LatticeFit_hexa(K,Wp,Wm,Wpperr,Wmmerr,Wpmerr,Echeck,x,l0,Weight)
[lp,lm]=LatticeEig_hexa(x,l0);

%fitting
d = 0.0;
Tot = sum(1-Echeck);
for i=1:length(K)
    if Echeck(i)==0
        if Weight==1 % standard weighted least square (more correlation than is currently accounted for?)
            D=((sqrt(lm(K(i,1),K(i,2)))-Wp(i))*Wmmerr(i)*(sqrt(lm(K(i,1),K(i,2)))-Wp(i))...
            -2*(sqrt(lm(K(i,1),K(i,2)))-Wp(i))*Wpmerr(i)*(sqrt(lp(K(i,1),K(i,2)))-Wm(i))...
            +(sqrt(lp(K(i,1),K(i,2)))-Wm(i))*Wpperr(i)*(sqrt(lp(K(i,1),K(i,2)))-Wm(i)))/(Wpperr(i)*Wmmerr(i)-Wpmerr(i)^2);
        else % unweighted
            D=((sqrt(lm(K(i,1),K(i,2)))-Wp(i)).^2+(sqrt(lp(K(i,1),K(i,2)))-Wm(i)).^2);
        end
        d=d+D/Tot/2;
        if D<0
            disp('Lattice_fit_hexa, failed');
            disp(x)
            disp([lm(K(i,1),K(i,2)),lp(K(i,1),K(i,2)),D/Tot/2])
        end
    end
end
d = sqrt(d); % rms error

end

