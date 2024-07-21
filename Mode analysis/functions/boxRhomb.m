function Box_rhomb = boxRhomb(b1, b2)
    BoxX=0.5.*[b1(1)+b2(1),b1(1)-b2(1),-b1(1)-b2(1),-b1(1)+b2(1),b1(1)+b2(1)];
    BoxY=0.5.*[b1(2)+b2(2),b1(2)-b2(2),-b1(2)-b2(2),-b1(2)+b2(2),b1(2)+b2(2)];

    Box_rhomb = [BoxX', BoxY'];
end