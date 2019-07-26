function [g,pos,h,p] = meshmesh(g,write,fname,fighnd,a,pos,hemisphere,affine,flip,inflate,checkori,dofillholes,optimise)
% Main brain-mesh plot. 
%
% Separates hemispheres, computes curvature, check orientation, inflate, 
% centre & fill holes. Plots the patch as a semi-transparent glass brain 
% for other plot functiopns to add to. 
%
% Returns everything in data.mesh
%

if iscell(pos)
   orig_pos = pos;
   pos      = pos{1};
end

if isempty(a);
    a = .6;
end

v = g.vertices;

% apply affine if req.
if length(affine) == 4
    fprintf('Applying affine transform\n');
    va = [v ones(length(v),1)]*affine;
    v  = va(:,1:3);
    g.vertices = v;
end

% flip x/y if required but POST affine transform
if flip
    v = g.vertices;
    v = v(:,[2 1 3]);
    g.vertices = v;
end


% check rotation
yl = max(v(:,2)) - min(v(:,2));
xl = max(v(:,1)) - min(v(:,1));

if xl > yl
    v       = v(:,[2 1 3]);
    g.faces = g.faces(:,[2 1 3]);
end

% centre and scale mesh
g.vertices = v - repmat(spherefit(v),[size(v,1),1]);


% ensure sourcemodel (pos) is around same scale as mesh boundaries
pos = fit_check_source2mesh(pos,g);

% explicitly optimise the alignment of the cloud points?
% (source model vertices and mesh vertices)
if optimise
    pos = align_clouds_3d(g.vertices,pos);
    clf; drawnow;
end

% m = min(g.vertices);% *1.1;
% M = max(g.vertices);% *1.1;
% 
% V        = pos - repmat(spherefit(pos),[size(pos,1),1]);
% V(:,1)   = m(1) + ((M(1)-m(1))).*(V(:,1) - min(V(:,1)))./(max(V(:,1)) - min(V(:,1)));
% V(:,2)   = m(2) + ((M(2)-m(2))).*(V(:,2) - min(V(:,2)))./(max(V(:,2)) - min(V(:,2)));
% V(:,3)   = m(3) + ((M(3)-m(3))).*(V(:,3) - min(V(:,3)))./(max(V(:,3)) - min(V(:,3)));
% pos      = V;

% % calculate curvature for shading
curv  = docurvature(struct('vertices',g.vertices,'faces',g.faces));
cvar  = var(curv);
dcurv = curv;

% inflate
if inflate
    fprintf('Inflating mesh\n'); nrep = 0;
    %while cvar / var(dcurv) < 50
        %nrep = nrep + 1;
        g = spm_mesh_inflate(struct('vertices',g.vertices,'faces',g.faces),200);
        %dcurv = docurvature(struct('vertices',g.vertices,'faces',g.faces));
    %end
    %fprintf('Finished inflation after %d reps\n',nrep);
end


if checkori
    b = CheckOrientationMesh(g);
    waitfor(b)
    g = b;
end


% only one hemisphere?
v = g.vertices;
f = g.faces;
c = spherefit(v);

% save centre point for later - e.g. label setting
g.centre = c;

left  = find(v(:,1) < c(1));
right = find(v(:,1) > c(1));

lfaces = find(sum(ismember(f,left),2)==3);
rfaces = find(sum(ismember(f,right),2)==3);

% return left/right indices
g.vleft            = v*NaN;
g.vleft(left,:)    = v(left,:);
g.vright           = v*NaN;
g.vright(right,:)  = v(right,:);
g.fleft            = f*NaN;
g.fleft(lfaces,:)  = f(lfaces,:);
g.fright           = f*NaN;
g.fright(rfaces,:) = f(rfaces,:);
g.curvature        = curv;



% fill holes stemming from hemisphere separation - if requested
%---------------------------------------------------------------
if dofillholes
    fprintf('Filling holes resulting from hemisphere separation..\n');
    
    % left side
    x = [];
    holecellarray = findTriMeshHoles(g.fleft,g.vleft);
    for i = 1:length(holecellarray)
        edges   = holecellarray{i};  % indices of vertices around edge
        ek      = nchoosek(edges,3); % triangles at all perimiter points
        x       = [x; ek];           % new faces
    end
    
    x0 = x;
    
    % update faces
    g.fleft = [g.fleft; x];
        
    % right side
    x = [];
    holecellarray = findTriMeshHoles(g.fright,g.vright);
    for i = 1:length(holecellarray)
        edges   = holecellarray{i};  % indices of vertices around edge
        ek      = nchoosek(edges,3); % triangles at all perimiter points
        x       = [x; ek];           % new faces
    end
    
    % update faces
    g.fright = [g.fright; x];
    
    % update full-model faces
    g.faces = [g.faces; x0; x];
end


switch hemisphere
    case{'left','L','l'}
        pg.vertices = g.vleft;
        pg.faces    = g.fleft;
        
        %pg.vertices         = v*NaN;
        %pg.vertices(left,:) = v(left,:);
        %pg.faces           = f*NaN;
        %pg.faces(lfaces,:) = f(lfaces,:);
  
    case{'right','R','r'}
        pg.vertices = g.vright;
        pg.faces    = g.fright;
        
        %pg.vertices          = v*NaN;
        %pg.vertices(right,:) = v(right,:);
        %pg.faces           = f*NaN;
        %pg.faces(rfaces,:) = f(rfaces,:);        
        
    otherwise
        pg = g;
end

% sanity check
if any( min(pg.faces(:)) == 0)
    bad = find(pg.faces == 0);
    pg.faces(bad) = NaN;
end

pg.faces = double(pg.faces);
pg.faces(pg.faces==0) = NaN;

% plot
if ~isempty(fighnd)
    if isnumeric(fighnd)
        % Old-type numeric axes handle
        h = plot(fighnd,gifti(pg));
    elseif ishandle(fighnd)
        % new for matlab2017b etc
        % [note editted gifti plot function]
        h = plot(gifti(pg),'fighnd',fighnd);
    end
else
    h  = plot(gifti(pg));
end
C = [.5 .5 .5];

set(h,'FaceColor',[C]); box off;
grid off;  set(h,'EdgeColor','none');
alpha(a); set(gca,'visible','off');

%h = get(gcf,'Children');
%set(h(end),'visible','off');
set(gca,'visible','off')

drawnow; hold on;

p = [];

if write == 1
    fprintf('Writing mesh gifti file: %s\n',[fname '.gii']);
    gout.vertices = g.vertices;
    gout.faces    = g.faces;
    gout = gifti(gout);
    save(gout,fname);
end



if exist('orig_pos')
    % is pos was a cell, re-instate 
    orig_pos{1} = pos;
    pos         = orig_pos;
end

end