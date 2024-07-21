function Eqs = modeEquations(U, Xref)
    % Equations
    Imag=0+1j;
    P = size(U);
    L2 = P(3);
    Num = P(1);

    uq=@(k1,k2,a) 1/sqrt(Num)*...
        transpose(exp(-1.0*Imag*(k1.*transpose(Xref(:,1))+...
                                 k2.*transpose(Xref(:,2))))*...
                  reshape(U(:,a,:),[P(1),P(3)]));
    M=@(k1,k2,a,b) sum(uq(-k1,-k2,a).*uq(k1,k2,b))/L2;
    SQ= @(k1,k2) sqrt(((M(k1,k2,1,1)+M(k1,k2,2,2))./2).^2-...
        (M(k1,k2,1,1).*M(k1,k2,2,2)-M(k1,k2,1,2).*M(k1,k2,2,1)));
    
    % Modes, 2 bands
    wp= @(k1,k2) (((M(k1,k2,1,1)+M(k1,k2,2,2))./2)+SQ(k1,k2)).^(-1/2);
    wm= @(k1,k2) (((M(k1,k2,1,1)+M(k1,k2,2,2))./2)-SQ(k1,k2)).^(-1/2);
        
    % Mode error, temporal part
    wpperr=@(k1,k2) sqrt(wp(k1,k2).^2 ./ 2 ./ (L2-1));
    wmmerr=@(k1,k2) sqrt(wm(k1,k2).^2 ./ 2 ./ (L2-1));
    wpmerr=@(k1,k2) sqrt(wm(k1,k2) .* 0);

    % Polarization
    yp= @(k1,k2) (((M(k1,k2,2,2)-M(k1,k2,1,1))./2)+SQ(k1,k2))./sqrt(M(k1,k2,1,2).*M(k1,k2,2,1)+(((M(k1,k2,2,2)-M(k1,k2,1,1))./2)+SQ(k1,k2)).^2);
    ym= @(k1,k2) (((M(k1,k2,2,2)-M(k1,k2,1,1))./2)-SQ(k1,k2))./sqrt(M(k1,k2,1,2).*M(k1,k2,2,1)+(((M(k1,k2,2,2)-M(k1,k2,1,1))./2)-SQ(k1,k2)).^2);
    xp= @(k1,k2) M(k1,k2,1,2)./sqrt(M(k1,k2,1,2).*M(k1,k2,2,1)+(((M(k1,k2,2,2)-M(k1,k2,1,1))./2)+SQ(k1,k2)).^2);
    xm= @(k1,k2) M(k1,k2,1,2)./sqrt(M(k1,k2,1,2).*M(k1,k2,2,1)+(((M(k1,k2,2,2)-M(k1,k2,1,1))./2)-SQ(k1,k2)).^2);
    
    %     Lp= @(k1,k2) (xp(k1,k2).*conj(xp(k1,k2)).*transpose(k1).^2+yp(k1,k2).^2.*transpose(k2).^2)./transpose(k1.^2+k2.^2);
    %     Lm= @(k1,k2) (xm(k1,k2).*conj(xm(k1,k2)).*transpose(k1).^2+ym(k1,k2).^2.*transpose(k2).^2)./transpose(k1.^2+k2.^2);
    Lm= @(k1,k2) (xm(k1,k2).*transpose(k1)+ym(k1,k2).*transpose(k2))...
              .*conj((xm(k1,k2).*transpose(k1)+ym(k1,k2).*transpose(k2)))...
              ./transpose(k1.^2+k2.^2);
    Lp= @(k1,k2) (xp(k1,k2).*transpose(k1)+yp(k1,k2).*transpose(k2))...
              .*conj((xp(k1,k2).*transpose(k1)+yp(k1,k2).*transpose(k2)))...
              ./transpose(k1.^2+k2.^2);

    Eqs.wp = wp;
    Eqs.wm = wm;
    Eqs.wpperr = wpperr;
    Eqs.wmmerr = wmmerr;
    Eqs.wpmerr = wpmerr;
    Eqs.Lm = Lm;
    Eqs.Lp = Lp;
    Eqs.M = M;
end