
% Example render AAL90 data on a pre-parcellated mesh
%-------------------------------------------------------------------

% Vector of values for each of the 90 regions in AAL90 order:
MyVector = randi([-50 50],90,1);

% render it using defaults: the quickest way
afigure; atemplate('overlay',MyVector,'method',{'aal_super',[]});