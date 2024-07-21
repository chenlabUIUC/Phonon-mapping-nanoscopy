function Path_hexa = pathHexa(b1, b2, SegN)
    % SegN = floor((n1 + n2)/2);
    
    Gamma = [0,0];
    M = [b2; b1+b2; b1]/2;
    R = zeros(4,2);
    R(1,:) = solveVertBisector(b2, -b1);
    R(2,:) = solveVertBisector(b2, b1+b2);
    R(3,:) = solveVertBisector(b1+b2, b1);
    R(4,:) = solveVertBisector(-b2, b1);
    KSegPaths = zeros(4,2,6);
    KSegPaths(:,:,1) = [Gamma; M(1,:); R(1,:); Gamma];
    KSegPaths(:,:,6) = [Gamma; M(1,:); R(2,:); Gamma];
    KSegPaths(:,:,3) = [Gamma; M(2,:); R(2,:); Gamma];
    KSegPaths(:,:,4) = [Gamma; M(2,:); R(3,:); Gamma];
    KSegPaths(:,:,5) = [Gamma; M(3,:); R(3,:); Gamma];
    KSegPaths(:,:,2) = [Gamma; M(3,:); R(4,:); Gamma];
    
    SegK_rec = [];
    SegDist_rec = [];
    
    for i = 1:6
        SegPt = KSegPaths(:,:,i);
        [SegK, SegDist] = KSegmentPath(SegPt, SegN);
        SegK_rec = cat(3,SegK_rec, SegK);
        SegDist_rec = cat(2,SegDist_rec, SegDist);
    end
    
    SegPtDist = mean(SegDist_rec(1:SegN-1:end,:),2);
    SegPtName = {'\Gamma','M','K','\Gamma'};

    Path_hexa.SegK_rec = SegK_rec;
    Path_hexa.SegDist_rec = SegDist_rec;
    Path_hexa.SegPtDist = SegPtDist;
    Path_hexa.SegPtName = SegPtName;
end


function p = solveVertBisector(v1,v2)
    M1 = v1/2;
    M2 = v2/2;
    slope1 = -1/(v1(2)/v1(1));
    slope2 = -1/(v2(2)/v2(1));
    syms x y
    eqn1 = y-M1(2) == slope1*(x-M1(1));
    eqn2 = y-M2(2) == slope2*(x-M2(1));
    [A,B] = equationsToMatrix([eqn1,eqn2],[x,y]);
    p = linsolve(A,B);
end

function [SegK, SegDist] = KSegmentPath(SegPt, SegN)
    for i = 1:size(SegPt,1)
        if i==1
            SegK = [0,0];
            SegDist = [0];
        else
            kx=linspace(real(SegPt(i-1,1)),real(SegPt(i,1)),SegN);
            ky=linspace(real(SegPt(i-1,2)),real(SegPt(i,2)),SegN);
            kxy = [kx;ky]';
            SegK = cat(1,SegK,kxy(2:end,:));
            Dist = sqrt(sum((kxy-SegPt(i-1,:)).^2,2))+SegDist(end);
            SegDist = cat(1,SegDist,Dist(2:end));
        end
    end
end


