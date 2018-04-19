function vol = sm2vol(sourcemodel,res,values,thresh)
%
%
%
%
% AS2018

if nargin < 2 || isempty(res);     res = [1 1 1]; end; if length(res) == 1; res = [res res res]; end
if nargin < 3 || isempty(values) ; dovals = 0;  else; dovals = 1; end
if nargin < 4  ; thresh = 10; else; thresh;     end % distance threshold


if ~isstruct(sourcemodel) 
    % OK if just passed vertex list, not fieldtrip sourcemodel structure
    v               = sourcemodel;
    sourcemodel     = struct;
    sourcemodel.pos = v;
    bounds = [floor(min(v)); ceil(max(v))];
    if res ~= 1
        sourcemodel.xgrid = linspace(bounds(1,1),bounds(2,1),res(1));
        sourcemodel.ygrid = linspace(bounds(1,2),bounds(2,2),res(2));
        sourcemodel.zgrid = linspace(bounds(1,3),bounds(2,3),res(3));
        sourcemodel.dim   = [length(sourcemodel.xgrid) length(sourcemodel.ygrid) length(sourcemodel.zgrid)];
    else
        sourcemodel.xgrid = bounds(1,1):bounds(2,1);
        sourcemodel.ygrid = bounds(1,2):bounds(2,2);
        sourcemodel.zgrid = bounds(1,3):bounds(2,3);
        sourcemodel.dim   = [length(sourcemodel.xgrid) length(sourcemodel.ygrid) length(sourcemodel.zgrid)];
    end
end
    
if res == 1 
    xgrid = sourcemodel.xgrid;
    ygrid = sourcemodel.ygrid;
    zgrid = sourcemodel.zgrid;
    dim   = sourcemodel.dim;
else
    xgrid = linspace(sourcemodel.xgrid(1),sourcemodel.xgrid(end),res(1));
    ygrid = linspace(sourcemodel.ygrid(1),sourcemodel.ygrid(end),res(2));
    zgrid = linspace(sourcemodel.zgrid(1),sourcemodel.zgrid(end),res(3));
    dim   = [length(xgrid) length(ygrid) length(zgrid)];
end


vol = zeros(dim);
pos = sourcemodel.pos;

numvox = 0;

% compute distance matrix
%dm     = cdist(pos,[xgrid' ygrid' zgrid']);
%n      = 0;

nit = prod(dim);
nit = round([0:.1:1]*nit);
nit = nit(2:end);

n     = 0;
for i = 1:dim(1)
    for j = 1:dim(2)
        for k = 1:dim(3)
            n = n + 1;
            if ismember(n,nit)
                fprintf('%d%% done\n',10*find(ismember(nit,n)));
            end
            
            xyz = [xgrid(i) ygrid(j) zgrid(k)];
            dm  = cdist(pos,xyz);
            
            [val,I] = min(dm);            
            if val  > thresh
                
            else
                % fill in this xyz location
                if dovals
                    this = values(I);
                    vol(i,j,k) = this;
                else
                    vol(i,j,k) = 1;
                end
                numvox = numvox + 1;
            end
            
        end
    end
end

fprintf('Volume has %d voxels\n',numvox);
