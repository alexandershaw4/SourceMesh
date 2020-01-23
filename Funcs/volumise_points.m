function vol = volumise_points(v,L)
%
% vol = volumise_points(v,[L])
%
% convert nx3 matrix, v, of 3D points to a 3D volume.
% optionally include functional value, L, for each point - must be length n
%
%
% AS 

dv = v;

if nargin < 2 || isempty(L)
    L = ones(length(v),1);
end

vol = zeros( (max(dv) - min(dv))+1 );
ndv = (min(dv))-1;

for i = 1:length(dv)
    if L(i) ~= 0
        a(1)  = L(i);
        a(2)  = vol(dv(i,1)-ndv(1),dv(i,2)-ndv(2),dv(i,3)-ndv(3));
        [~,I] = max(abs(a));
        vol(dv(i,1)-ndv(1),dv(i,2)-ndv(2),dv(i,3)-ndv(3)) = ...
            vol(dv(i,1)-ndv(1),dv(i,2)-ndv(2),dv(i,3)-ndv(3)) + a(I);
    end
end