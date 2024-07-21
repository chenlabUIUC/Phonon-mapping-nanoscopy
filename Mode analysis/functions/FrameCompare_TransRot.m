function C=FrameCompare_TransRot(x,Q)
    % Loss function for rotation correction
    % Q is particle positions by frame (n,d,t)
    % x is array of rotation angle of each frame (t,)
    
    %x[1] is assumed/forced to be zero
    x(1,:)=0;
    C=0;

    rot_mat = @(theta) [cos(theta) sin(theta);
                       -sin(theta) cos(theta)];

    % sum of mean squared distance between all pairs of frames.
    t = size(Q,3);

    for i = 1:t
        pos1 = Q(:,:,i) - x(i,2:3);
        for j = i+1:t
            pos2 = Q(:,:,j) - x(j,2:3);
            C = C + 2*sum((pos1*rot_mat(x(i,1))-pos2*rot_mat(x(j,1))).^2,'all')/t.^2;
        end
    end

end