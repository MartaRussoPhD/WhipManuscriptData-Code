function colNew = ColorModify(colSource,H,S,L)
%Modify color by HSL
    tempCol = rgb2hsv(colSource);
    tempCol = tempCol .* [H,S,L];
    tempCol(tempCol>1) = 1;
    colNew = hsv2rgb(tempCol);
end