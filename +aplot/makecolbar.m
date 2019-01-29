function RGB = makecolbar(I,netcmap)
% Register colorbar values to an overlay /  T-vector
%

if any(any(netcmap ~= 0))
    Colors = colormap(netcmap);
else
    Colors   = jet;
end

NoColors = length(Colors);

Ireduced = (I-min(I))/(max(I)-min(I))*(NoColors-1)+1;
RGB      = interp1(1:NoColors,Colors,Ireduced);

end