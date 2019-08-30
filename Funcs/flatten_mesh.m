function [dv,uv,g,bnd] = flatten_mesh(v,f,c,doplot)
% squash a 3D cloud [input v(nx3) ] into a 2D orthographic projection, 
% using x/y spread relative to original depth. for brains.
%
%  % squash/flatten v(nx3) to dv(nx2):
%  [dv,uv,g] = flatten_mesh(v)
%
%  % approximately translate back to 3D:
%  ov      = dv./uv;     % first 2 dims in orig space
%  ov(:,3) = diff(uv')'; % approx (quantised) 3rd dim [scale is wrong]
%
% AS

if nargin < 4
    doplot=0;
end

% estimate size of required 2d canvas
v      = v - repmat(spherefit(v),[size(v,1),1]);
bnd    = [min(v); max(v)];
d3s    = diff(bnd(:,3)) / 2;
newcnv = [bnd(1,1:2) - d3s ; bnd(2,1:2) + d3s];
former = bnd(1:2,1:2);

% sort highest to lowest points
[~,top2b] = sort(v(:,3),'descend');

n_g = 5000; % bigger=smoother

% chunk top-to-bottom vertex list into n_g groups, 1:n_g = incrs. depth
depths = linspace(1,length(top2b),n_g+1);
dgroup = zeros(1,length(top2b));
for i = 1:n_g
    dgroup(depths(i):depths(i+1)) = i;
end

% depth group indexing of the original vertex sequence:
g(top2b)=dgroup;

% scale factor determines how far we spread in x and y directions as a
% function of depth (former z direction)
scale_factor = ( newcnv./former ) / n_g;
scale_factor = mean(scale_factor,1);
new_list     = [];

for i   = 1:length(v)
    pnt = v(i,:);
    grp = g(i);
    
    new_pnt = pnt(1:2) .* (grp*scale_factor);
    new_list = [new_list ; new_pnt];
    
    % retain multiplier (e.g. dv = v(1:2).*uv]
    uv(i,:) = (grp*scale_factor);
end

% complete 2nd vertex list - with vertices in same order as input
dv = new_list;

if doplot
    % plot example - brain data:
    %[dv,uv] = flatten_mesh(v);
    if ~isempty(f)
        p = struct('faces',[f(:,1:2)],'vertices',[dv]);
        p = patch(p);
        set(p,'EdgeColor',[.4 .4 .4]);p.EdgeAlpha=.4; hold on;
    end
    
    sc = scatter(dv(:,1),dv(:,2),80,c,'filled');hold on
    sc.MarkerFaceAlpha = .4;
    k  = boundary(double(v(:,1:2)));
    k0 = boundary(double(dv));
    l0 = line(v(k,1),v(k,2));
    l1 = line(dv(k0,1),dv(k0,2));
    
    l0.Color = [.4 .4 .4];
    l0.LineWidth = 2;
    l1.Color = [.4 .4 .4];
    l1.LineWidth = 2;
    axis square;
end



