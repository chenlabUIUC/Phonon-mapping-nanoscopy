function standardModePlot(K, W, cmp, xy_range, c_range, pointSize, title_text)
    scatter3(K(:,1),K(:,2),W, pointSize,W,'filled'); 
    view(2)
    daspect([1 1 1]); box on; axis equal; hold on;
    colormap(gca, cmp)
    clim([0, c_range])
    axis([-1,1,-1,1] .* xy_range);
    title(title_text)
end