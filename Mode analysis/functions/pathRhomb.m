function Path_rhomb = pathRhomb(b1, b2, SegN)
    % SegN = floor((n1 + n2)/2);

    Gamma = [0,0]; M = b1/2; R1 = (b1+b2)/2; R2 = (b1-b2)/2;
    R1mag = sqrt(sum(R1.^2)); R2mag = sqrt(sum(R2.^2));
    if R1mag>R2mag; R3=R1; R1=R2; R2=R3; end
    R1mag = sqrt(sum(R1.^2)); R2mag = sqrt(sum(R2.^2));
    Mmag = sqrt(sum(M.^2)); R1Mmag = sqrt(sum((R1-M).^2)); R2Mmag = sqrt(sum((R2-M).^2));

    normVec = [R1(2),-R1(1)]; normVecMag = sqrt(sum(normVec.^2));
    normVec = normVec/normVecMag;
    ReflectMat = [1-2*normVec(1)^2,-2*normVec(1)*normVec(2);...
                  -2*normVec(1)*normVec(2), 1-2*normVec(2)^2];

    SegPt = [Gamma;R1;M;Gamma;R2;M]; % Main.
    KSegPaths = zeros(6,2,2);
    KSegPaths(:,:,1) = [Gamma;R1;M;Gamma;R2;M];
    KSegPaths(:,:,2) = [Gamma;R1;M;Gamma;R2;M] * ReflectMat;

    SegK_rec = [];
    SegDist_rec = [];
    
    for i = 1:size(KSegPaths,3)
        SegPt = KSegPaths(:,:,i);
        [SegK, SegDist] = KSegmentPath(SegPt, SegN);
        SegK_rec = cat(3,SegK_rec, SegK);
        SegDist_rec = cat(2,SegDist_rec, SegDist);
    end

    SegPtDist = SegDist(1:SegN-1:end);
    SegPtName = {'\Gamma','R_1','M','\Gamma','R_2','M'};

    Path_rhomb.SegK_rec = SegK_rec;
    Path_rhomb.SegDist_rec = SegDist_rec;
    Path_rhomb.SegPtDist = SegPtDist;
    Path_rhomb.SegPtName = SegPtName;
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


