function Box_hexa = boxHexa(b1, b2)
    BoxN_hexa = [1/3,-1/3; 2/3,1/3; 1/3,2/3;...
                 -1/3,1/3; -2/3,-1/3; -1/3,-2/3];
    Box_hexa = b1'*BoxN_hexa(:,1)' + b2'*BoxN_hexa(:,2)';
    Box_hexa = Box_hexa';
    Box_hexa = Box_hexa([1:end,1],:);
end