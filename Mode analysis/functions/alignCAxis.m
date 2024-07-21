function alignCAxis(n)
    cmin = 0;
    cmax = -inf;
    for i = n
        subplot(3,4,i)
        clim auto
        c = clim();
        if c(2)>cmax
            cmax = c(2);
        end
    end
    for i = n
        subplot(3,4,i)
        clim([cmin, cmax])
    end
end