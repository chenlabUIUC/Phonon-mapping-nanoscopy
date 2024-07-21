function PlotPositionMatrix(pos)
    temp = permute(pos,[1,3,2]);
    plot(temp(:,:,1)', temp(:,:,2)','o-');
end