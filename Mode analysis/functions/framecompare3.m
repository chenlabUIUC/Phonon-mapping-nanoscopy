function C=framecompare2(x,Q)

    % Q is particle positions by frame (n,d,t)
    % x is array of rotation angle of each frame (t,)
    
    %x[1] is assumed/forced to be zero
    x(1)=0;
    C=0;

    rot_mat = @(theta) [cos(theta) sin(theta);
                       -sin(theta) cos(theta)];

    % sum of mean squared distance between any 2 frames.
    for i=1:size(Q,3)
        pos1 = reshape(Q(:,:,i),size(Q,1),[]);
        for j=i+1:size(Q,3)
            pos2 = reshape(Q(:,:,j),size(Q,1),[]);
            C = C + 2*sum(sum((pos1*rot_mat(x(i))-pos2*rot_mat(x(j))).^2))/size(Q,3).^2;
        end
    end

end