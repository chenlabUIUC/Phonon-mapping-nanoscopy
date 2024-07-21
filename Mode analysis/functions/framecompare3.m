function C = framecompare3(x,Q)

    % Q is particle positions by frame (n,d,t)
    % x is array of rotation angle of each frame (t,)
    
    %x[1] is assumed/forced to be zero
    x(1)=0;
    C=0;

    rot_mat = @(theta) [cos(theta) sin(theta);
                       -sin(theta) cos(theta)];

    for i = 1:size(Q,3)
        Q(:,:,i) = Q(:,:,i) * rot_mat(x(i));
    end

    t = size(Q,3);

    % sum of mean squared distance between any 2 frames.
    for i=1:t
        pos1 = Q(:,:,i);
        temp = Q - pos1;
        C = C + sum(temp.^2,'all')/t/t;
    end

end