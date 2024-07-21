function [b1, b2] = reciprocalLattice(a1, a2)
    % Extract components of a1 and a2
    a1x = a1(1);
    a1y = a1(2);
    a2x = a2(1);
    a2y = a2(2);
    
    % Calculate the area of the unit cell in real space
    A = a1x * a2y - a1y * a2x;
    
    % Check for zero area (linearly dependent vectors)
    if A == 0
        error('The lattice vectors a1 and a2 are linearly dependent or the area is zero.');
    end
    
    % Calculate the reciprocal lattice vectors
    b1 = (2 * pi / A) * [a2y, -a2x];
    b2 = (2 * pi / A) * [-a1y, a1x];
end