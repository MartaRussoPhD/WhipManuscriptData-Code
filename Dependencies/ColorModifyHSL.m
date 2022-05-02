function colNew = ColorModifyHSL(colSource,H,S,L)
    tempCol = rgb2hsv(colSource);
    tempCol = tempCol .* [H,S,L];
    tempCol(tempCol>1) = 1;
    colNew = hsv2rgb(tempCol);
end