function y = NewMeanFilt3D(M,m1,m2,m3)
% Mean window filtering for 3D matrices
%
% 
%
% AS2016 [updt] [updt 17]


if ndims(M) > 3; fprintf('Too many dimensions!\n'); return; end
if ndims(M) < 3; fprintf('Too few dimensions!\n'); return; end

if nargin < 4
    [m2, m3] = deal(m1);
end

y = smoothmat3(M,m1,m2,m3);

end

function matrixOut = smoothmat3(matrixIn,Nr,Nc,Nz)

N(1) = Nr;  N(2) = Nc;  N(3) = Nz;

[row,col,dep] = size(matrixIn);

eL = spdiags(ones(row,2*N(1)+1),(-N(1):N(1)),row,row);
eR = spdiags(ones(col,2*N(2)+1),(-N(2):N(2)),col,col);
eZ = spdiags(ones(dep,2*N(3)+1),(-N(3):N(3)),dep,dep);

A = isnan(matrixIn);
matrixIn(A) = 0;

for i = 1:length(eZ);        
    nrm(:,:,i) = eL*(~A(:,:,i))*eR;
    nrm(:,:,i) = nrm(:,:,i)*eZ(i,i);
end

nrmlize    = nrm;
nrmlize(A) = NaN;

for i = 1:length(eZ);
    matrixOut(:,:,i) = eL*matrixIn(:,:,i)*eR;
    matrixOut(:,:,i) = matrixOut(:,:,i)*eZ(i,i);
end
matrixOut = matrixOut./nrmlize;

end

