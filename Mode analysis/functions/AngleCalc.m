function ang = AngleCalc(p0,p1,p2)
    if nargin<3
        temp = zeros(size(p0));
        temp(1) = 1;
        p2 = p0 + temp;
    end

    n1 = (p2 - p0) / norm(p2 - p0);  % Normalized vectors
    n2 = (p1 - p0) / norm(p1 - p0);

    ang = atan2(det([n2; n1]), dot(n1, n2));
    % ang = acos(dot(n1,n2)/norm(n1)/norm(n2));
    ang = ang/pi*180;
end