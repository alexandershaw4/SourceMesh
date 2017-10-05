function varargout = atemplate(varargin)
% Add AAL90 networks and overlays to a smoothed brain mesh
%
% NOTE: Unknown BUG when using with Matlab 2015a on Linux.
% Working with Matlab 2014a & 2017a on Mac & Matlab 2012a on Linux
%
% Does this:
%
% 'overlay', plot functional overlay
% 'network', plot nodes & edges from AAL network analysis
% 'nodes'  , plot only nodes from AAL90 atlas
% 'tracks' , tracks, header (as read with trk_read^ )
% 'labels' , add AAL90 labels
% 'video'  , make video overlay 
% 
%
%
%  Usages:
%
%  atemplate('labels');         template mesh with AAL labels
%  atemplate('overlay',L);      template mesh with overlay [90x1]
%  atemplate('network',A);      template mesh with nodes & edges [90x90]
%  atemplate('nodes', N);       template mesh with nodes labelled [90,1]
%  atemplate('overlay',L,'network',A,'labels'); overlay, network & labels
%  atemplate('tracks',tracks,header); plot tracks loaded with trk_read
%  atemplate('gifti',g);        use supplied gifti surface / mesh 
%
%  atemplate('gifti'  ,g,'write',name);  write mesh gifti
%  atemplate('overlay',L,'write',name);  write mesh & overlay giftis
%
%  atemplate('gifti',g,'video',m,'name',times); m = [90 * ntimes] 
%                                               name  = video savename
%                                               times = vector of titles
%
%  Note: any combination of the inputs is possible.
%
%
%  Cortical mesh from mri
%  ---------------------------
%  If the gifti option is included, input (g) may be the filename of a 
%  coregistered ctf .mri file. This will call Vol2SurfAS which uses 
%  fieldtrip & isosurface to normalise, align, segment and extract a
%  cortical surface. This is then centred, smoothed and converted to a
%  gifti object.
%
%  See also: slice3() slice2()
%
% ^trk_read requires along-tract-stats toolbox
%
% AS17

pmesh  = 1;
labels = 0;
write  = 0;
fname  = [];
fighnd = [];
colbar = 1;
for i  = 1:length(varargin)
    if strcmp(varargin{i},'overlay'); L = varargin{i+1};     end
    if strcmp(varargin{i},'network'); A = varargin{i+1};     end
    if strcmp(varargin{i},'tracks');  T = varargin{i+1}; H = varargin{i+2}; end
    if strcmp(varargin{i},'labels');  labels = 1;            end
    if strcmp(varargin{i},'nosurf');  pmesh  = 0;            end
    if strcmp(varargin{i},'nodes');   N = varargin{i+1};     end
    if strcmp(varargin{i},'gifti');   g = varargin{i+1};     end
    if strcmp(varargin{i},'write');   write  = 1; fname = varargin{i+1}; end
    if strcmp(varargin{i},'fighnd');  fighnd = varargin{i+1}; end
    if strcmp(varargin{i},'nocolbar');colbar = 0; end
    if strcmp(varargin{i},'video');   V     = varargin{i+1};...
                                      fpath = varargin{i+2};
                                      times = varargin{i+3}; end
    if strcmp(varargin{i},'othermesh'); M = varargin{i+1};     
                                        O = varargin{i+2};   end    
    
end


% Plot Surface
%-------------------------------------------------------------
try   mesh = g;
      fprintf('User provided mesh\n');
      if ischar(mesh);
          % generate a cortical mesh from the mri (Vol2SurfAS)
          fprintf('MRI (%s) is character (a filename?): attempting to load, segment & mesh\n',mesh);
          mesh = Vol2SurfAS(mesh,'ctf','smooth',0.15);
          fprintf('\n\nSuccessfully generated subject mesh from .mri!\n');
      end
catch mesh = read_nv();
      fprintf('Template mesh\n');
end


% plot the glass brain
if     pmesh && ~exist('T','var');
       mesh = meshmesh(mesh,write,fname,fighnd);
elseif pmesh
       mesh = meshmesh(mesh,write,fname,fighnd,.3);
end


% find closest vertices and overlay
try L; overlay(mesh,L,write,fname,colbar);end 

% o colbar is plotting both overlay & network!
isover = exist('L','var') || exist('V','var');
if  isover && exist('A','var') 
    colbar = 0;
    alpha(.2);
end


% draw edges and edge-connected nodes
try A; connections(A,colbar);      end 
% draw dti tracks loaded with trk_read
try T; drawtracks(T,H,mesh);       end 
% draw N(i) = 1 nodes
try N; drawnodes(N);               end 


try M;
    [M,y] = othermesh(M,O,mesh,write,fname);
    varargout{1} = M;
    varargout{2} = y;
end

if labels; 
    if     exist('A','var'); addlabels(A);
    elseif exist('N','var'); addlabels(diag(N));
    else   addlabels(ones(90,90));
    end
end

try V; 
    tv = 1:size(V,2);
    try tv = times; end
    video(mesh,V,1,fpath,tv); 
end


end



function connections(A,colbar)
% Network (Node & Edges) plotter.
%
%


% Edges
%-------------------------------------------------------------
[node1,node2,strng] = conmat2nodes(A);
RGB = makecolbar(strng);

% LineWidth (scaled) for strength
if any(strng)
    R = [min(strng),max(strng)];
    S = ( strng - R(1) ) + 1e-3;
    
    % If all edges same value, make thicker
    if  max(S(:)) == 1e-3; 
        S = 3*ones(size(S)); 
    end
else
    S = [0 0];
end

% If too few strengths, just use red edges
LimC = 1;
if all(all(isnan(RGB)))
    RGB  = repmat([1 0 0],[size(RGB,1) 1]);
    LimC = 0;
end

% Paint edges

for i = 1:size(node1,1)
    line([node1(i,1),node2(i,1)],...
        [node1(i,2),node2(i,2)],...
        [node1(i,3),node2(i,3)],...
        'LineWidth',S(i),'Color',[RGB(i,:)]);
end

% Set colorbar only if there are valid edges
if any(i) && colbar
    set(gcf,'DefaultAxesColorOrder',RGB)
    if colbar
        colorbar
    end
end
if LimC && colbar
    caxis(R);
end

drawnow;


% Nodes (of edges only)
%-------------------------------------------------------------
hold on;
for i = 1:size(node1,1)
    scatter3(node1(i,1),node1(i,2),node1(i,3),'filled','k');
    scatter3(node2(i,1),node2(i,2),node2(i,3),'filled','k');
end

drawnow;

end

function drawnodes(N)
% Node plotter. N = (90,1) with 1s for nodes to plot and 0s to ignore.
%
% 

hold on;
% load([fileparts(which('conmat2nodes')),'/AAL_SOURCEMOD.mat']);
% v = template_sourcemodel.pos;
% align aal labels to mesh
v = fixmesh(mesh); % get AAL positions on mesh surface


ForPlot = v(find(N),:);

for i = 1:length(ForPlot)
    scatter3(ForPlot(i,1),ForPlot(i,2),ForPlot(i,3),'filled','r');
end


end

function RGB = makecolbar(I)
% Register colorbar values to our overlay /  T-vector
%

Colors   = jet;
NoColors = length(Colors);

Ireduced = (I-min(I))/(max(I)-min(I))*(NoColors-1)+1;
RGB      = interp1(1:NoColors,Colors,Ireduced);

end

function drawtracks(tracks,header,mesh)
% IN PROGRESS - NOT SURE HOW TO DEAL WITH TRACKS
%
% - Use trk_read from 'along-tract-stats' toolbox
%

hold on; clc;
All = [];

% put all tracks into a single matrix so we can fit a sphere
for iTrk = 1:length(tracks)
    if iTrk > 1; fprintf(repmat('\b',size(str))); end
    str = sprintf('Building volume for sphere fit (%d of %d)\n',iTrk,length(tracks));
    fprintf(str);
    
    matrix = tracks(iTrk).matrix;
    matrix(any(isnan(matrix(:,1:3)),2),:) = [];
    All = [All ; matrix];
end

% centre on 0 by subtracting sphere centre
iAll      = All;
iAll(:,1) = All(:,1)*-1;
iAll(:,2) = All(:,2)*-1;
Centre = spherefit(iAll);
maxpts = max(arrayfun(@(x) size(x.matrix, 1), tracks));

% Use minmax template vertices as bounds
MM(1,:) = min(mesh.vertices);
MM(2,:) = max(mesh.vertices);
MT(1,:) = min(iAll-repmat(Centre,[size(iAll,1),1]));
MT(2,:) = max(iAll-repmat(Centre,[size(iAll,1),1]));

pullback = min(MM(:,2)) - min(MT(:,2));
pullup   = max(MM(:,3)) - max(MT(:,3));

D = mean((MM)./MT);

% this time draw the tracks
for iTrk = 1:length(tracks)
    matrix = tracks(iTrk).matrix;
    matrix(any(isnan(matrix(:,1:3)),2),:) = [];
    
    matrix(:,1) = matrix(:,1)*-1; % flip L-R
    matrix(:,2) = matrix(:,2)*-1; % flip F-B
    M           = matrix - repmat(Centre,[size(matrix,1),1]); % centre
    M           = M.*repmat(D,[size(M,1),1]);
    M(:,2)      = M(:,2) + (pullback*1.1);                          % pullback
    M(:,3)      = M(:,3) + (pullup*1.1);                            % pull up

    h = patch([M(:,1)' nan], [M(:,2)' nan], [M(:,3)' nan], 0);
    cdata = [(0:(size(matrix, 1)-1))/(maxpts) nan];
    set(h,'cdata', cdata, 'edgecolor','interp','facecolor','none');
end

h = get(gcf,'Children');
set(h,'visible','off');

end


function overlay(mesh,L,write,fname,colbar)
% Functional overlay plotter
%
% mesh is the gifti / patch
% L is the overlay (90,1)
% write is boolean flag
% fname is filename is write = 1;
%


% interp shading between nodes or just use mean value?
interpl = 1; 

% Overlay
%load('AAL_SOURCEMOD');          % get AAL source vertices
%v  = template_sourcemodel.pos;  % AAL vertices
v = fixmesh(mesh); % get AAL positions on mesh surface

x  = v(:,1);                    % AAL x verts
mv = mesh.vertices;             % brain mesh vertices
nv = length(mv);                % number of brain vertices
OL = sparse(length(L),nv);      % this will be overlay matrix we average
r  = 1200;                      % radius - number of closest points on mesh
w  = linspace(.1,1,r);          % weights for closest points
w  = fliplr(w);                 % 
M  = zeros( length(x), nv);     % weights matrix: size(len(mesh),len(AAL))


% if is same verts as mri, just rescale & overlay
S  = [min(L(:)),max(L(:))];
% otherwise find closest points (assume both in mm)
fprintf('Determining closest points between AAL & template vertices (overlay)\n');
for i = 1:length(x)
    
    % reporting
    if i > 1; fprintf(repmat('\b',[size(str)])); end
    str = sprintf('%d/%d',i,(length(x)));
    fprintf(str);
    
    % find closest point[s] in cortical mesh
    dist       = cdist(mv,v(i,:));
    [junk,ind] = maxpoints(dist,r,'min');
    OL(i,ind)  = w*L(i);
    M (i,ind)  = w;
    
end
fprintf('\n');

if ~interpl
     % mean value of a given vertex
    OL = mean((OL),1);
else
    for i = 1:size(OL,2)
        % average overlapping voxels
        L(i) = sum( OL(:,i) ) / length(find(OL(:,i))) ;
    end
    OL = L;
end


% normalise and rescale
y  = S(1) + ((S(2)-S(1))).*(OL - min(OL))./(max(OL) - min(OL));

y(isnan(y)) = 0;
y  = full(y);

% spm mesh smoothing
fprintf('Smoothing overlay...\n');
y = spm_mesh_smooth(mesh, y(:), 4);
hh = get(gca,'children');
set(hh(end),'FaceVertexCData',y(:),'FaceColor','interp');
drawnow;
shading interp
colormap('jet');

if colbar
    drawnow; pause(.5);
    colorbar('peer',gca,'South');    
end
    
    if write;
        fprintf('Writing overlay gifti file: %s\n',[fname 'Overlay.gii']);
        g       = gifti;
        g.cdata = double(y);
        g.private.metadata(1).name  = 'SurfaceID';
        g.private.metadata(1).value = [fname 'Overlay.gii'];
        save(g, [fname  'Overlay.gii']);
    end
    
end



% function [M,y] = othermesh(M,O,mesh,write,fname)
% 
% colbar = 1;
% 
% fprintf('Fitting supplied mesh vertices to template mesh\n');
% nv = mesh.vertices;
% ov = M.vertices;
% 
% %nv = nv - repmat(spherefit(nv),[size(nv,1),1]);
% ov = ov - repmat(spherefit(ov),[size(ov,1),1]);
% 
% bounds = [min(nv); max(nv)];
% orig   = [min(ov); max(ov)];
% m      = mean( bounds./orig );
% new    = (0.99)* ( ov.*repmat(m,[size(ov,1),1]) );
% b      = bounds;
% 
% new(:,1) = b(1,1) + ((b(2,1)-b(1,1))) .* (new(:,1) - min(new(:,1)))./(max(new(:,1))-min(new(:,1)));
% new(:,2) = b(1,2) + ((b(2,2)-b(1,2))) .* (new(:,2) - min(new(:,2)))./(max(new(:,2))-min(new(:,2)));
% new(:,3) = b(1,3) + ((b(2,3)-b(1,3))) .* (new(:,3) - min(new(:,3)))./(max(new(:,3))-min(new(:,3)));
% 
% 
% L  = O.cdata;
% x  = new(:,1);
% mv = mesh.vertices;
% v  = new;
% OL = sparse(length(L),length(mv)); % this will be overlay matrix we average
% 
% r  = ceil(length(mv)*0.0705);      % radius - number of closest points on mesh
% 
% w  = linspace(.1,1,r);          % weights for closest points
% w  = fliplr(w);                 % 
% M  = zeros( length(x), length(mesh.vertices) );
% %w  = ones(size(w)); 
% 
% S  = [min(L(:)),max(L(:))];
% % otherwise find closest points (assume both in mm)
% fprintf('Determining closest points between supplied mesh & template vertices (overlay)\n');
% for i = 1:length(x)
%     
%     % reporting
%     if i > 1; fprintf(repmat('\b',[size(str)])); end
%     str = sprintf('%d/%d',i,(length(x)));
%     fprintf(str);
%     
%     % find closest point[s] in cortical mesh
%     dist       = cdist(mv,v(i,:));
%     [junk,ind] = maxpoints(dist,r,'min');
%     OL(i,ind)  = w*L(i);
%     M (i,ind)  = w;
%     
% end
% fprintf('\n');
% 
% OL = killinterhems(OL);
% M  = killinterhems(M);
% 
% interpl = 1;
% L = double(L);
% if ~interpl
%     OL = mean((OL),1); % mean value of a given vertex
% else
%     for i = 1:size(OL,2)
%         % average overlapping voxels
%         L(i) = sum( OL(:,i) ) / length(find(OL(:,i))) ;
%     end
%     OL = L;
%     %OL = full(L);
% end
% 
% % normalise and rescale
% y  = S(1) + ((S(2)-S(1))).*(OL - min(OL))./(max(OL) - min(OL));
% %y  = S(1) + ((S(2)-S(1))).*TSNorm(OL); % or normalise as sparse
% 
% y(isnan(y)) = 0;
% y  = full(y);
% y = double(y);
% 
% 
% % spm mesh smoothing
% fprintf('Smoothing overlay...\n');
% y = spm_mesh_smooth(mesh, y(:), 4);
% 
% y  = y(:);
% hh = get(gca,'children');
% 
% set(hh(end),'FaceVertexCData',y, 'FaceColor','interp');
% shading interp
% if colbar
%     colorbar
% end
% 
% if write;
%     fprintf('Writing overlay gifti file: %s\n',[fname 'Overlay.gii']);
%     g       = gifti;
%     g.cdata = double(y);
%     g.private.metadata(1).name  = 'SurfaceID';
%     g.private.metadata(1).value = [fname 'Overlay.gii'];
%     save(g, [fname  'Overlay.gii']);
% end
% 
% 
% 
% 
% end

function x = killinterhems(x);

S  = size(x);
xb = (S(1)/2)+1:S(1);
yb = (S(2)/2)+1:S(2);
xa = 1:S(1)/2;
ya = 1:S(2)/2;

x(xa,yb) = 0;
x(xb,ya) = 0;

end




function newpos = fixmesh(g)
% plot as transparent grey gifti surface
%
% AS

v = g.vertices;
v = v - repmat(spherefit(v),[size(v,1),1]);% Centre on ~0
g.vertices=v;


%
load('AAL_SOURCEMOD');
pos  = template_sourcemodel.pos;


% Centre on ~0
pos = pos - repmat(spherefit(pos),[size(pos,1),1]);

for i = 1:length(pos)
    
    this = pos(i,:);
    [t,I] = maxpoints(cdist(v,this),1,'max');
    
    newpos(i,:) = v(I,:);
    
end



% M        = max(v);
% m        = min(v);
% 
% M = M*1.15;
% m = m*1.15;
% 
% V        = g.vertices;
% V        = V - repmat(spherefit(V),[size(V,1),1]);
% 
% V(:,1)   = m(1) + ((M(1)-m(1))).*(V(:,1) - min(V(:,1)))./(max(V(:,1)) - min(V(:,1)));
% V(:,2)   = m(2) + ((M(2)-m(2))).*(V(:,2) - min(V(:,2)))./(max(V(:,2)) - min(V(:,2)));
% V(:,3)   = m(3) + ((M(3)-m(3))).*(V(:,3) - min(V(:,3)))./(max(V(:,3)) - min(V(:,3)));
% 
% g.vertices = V;

end

function g = meshmesh(g,write,fname,fighnd,a);

if nargin < 5;
    a = .6;
end

% plot
if ~isempty(fighnd)
    %axes(fighnd);
    h = plot(fighnd,gifti(g));
else
    h = plot(gifti(g));
end
C = [.5 .5 .5];

set(h,'FaceColor',[C]); box off;
grid off;  set(h,'EdgeColor','none');
alpha(a); set(gca,'visible','off');

h = get(gcf,'Children');
set(h(end),'visible','off');
drawnow;

if write;
    fprintf('Writing mesh gifti file: %s\n',[fname '.gii']);
    g = gifti(g);
    save(g,fname);
end


end

function addlabels(V)
% Add AAL labels to the plot
%
% Finds indices of V~=0, which is (90,1) || (1,90)
%
%

load('labels');
labels = strrep(labels,'_',' ');

%load('AAL_SOURCEMOD');
pos = fixmesh(mesh); % get AAL positions on mesh surface


% compile list of in-use node indices
%------------------------------------
to = []; from = [];
for i  = 1:size(V,1)
    ni = find(logical(V(i,:)));
    if any(ni)
        to   = [to   ni];
        from = [from repmat(i,[1,length(ni)]) ];
    end
end

AN  = unique([to,from]);
%v   = template_sourcemodel.pos;
off = 1.5;

% add these to plot with offset
%------------------------------------
for i = 1:length(AN)
    L = labels{AN(i)};
    switch L(end)
        case 'L';
            t(i) = text(v(AN(i),1)-(off*5),v(AN(i),2)-(off*5),v(AN(i),3)+off,L);
        case 'R';
            t(i) = text(v(AN(i),1)+(off*2),+v(AN(i),2)+(off*2),v(AN(i),3)+off,L);
    end
end
set(t,'Fontsize',10)

end

function Centre = spherefit(X)
% Fit sphere to centre of vertices, return centre points
%
%

A =  [mean(X(:,1).*(X(:,1)-mean(X(:,1)))), ...
    2*mean(X(:,1).*(X(:,2)-mean(X(:,2)))), ...
    2*mean(X(:,1).*(X(:,3)-mean(X(:,3)))); ...
    0, ...
    mean(X(:,2).*(X(:,2)-mean(X(:,2)))), ...
    2*mean(X(:,2).*(X(:,3)-mean(X(:,3)))); ...
    0, ...
    0, ...
    mean(X(:,3).*(X(:,3)-mean(X(:,3))))];
A = A+A.';
B = [mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,1)-mean(X(:,1))));...
     mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,2)-mean(X(:,2))));...
     mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,3)-mean(X(:,3))))];
Centre=(A\B).';
end


function video(mesh,L,colbar,fpath,tv)
%

num = 1; % number of brains

% interp shading between nodes or just use mean value?
interpl = 1; 

% Overlay
%load('AAL_SOURCEMOD');          % get AAL source vertices
%v  = template_sourcemodel.pos;  % AAL vertices
v = fixmesh(mesh); % get AAL positions on mesh surface

x  = v(:,1);                    % AAL x verts
mv = mesh.vertices;             % brain mesh vertices
nv = length(mv);                % number of brain vertices
ntime = size(L,2);
OL = zeros(size(L,1),nv,ntime);      % this will be overlay matrix we average
r  = 1200;                      % radius - number of closest points on mesh
w  = linspace(.1,1,r);          % weights for closest points
w  = fliplr(w);                 % 
M  = zeros( length(x), nv);     % weights matrix: size(len(mesh),len(AAL))


S  = [min(L)',max(L)'];

% find closest points (assume both in mm)
fprintf('Determining closest points between AAL & template vertices (overlay)\n');
for i = 1:length(x)

    % reporting
    if i > 1; fprintf(repmat('\b',[size(str)])); end
    str = sprintf('%d/%d',i,(length(x)));
    fprintf(str);    

    % find closest point[s] in cortical mesh
    dist       = cdist(mv,v(i,:));
    [junk,ind] = maxpoints(dist,r,'min');
    OL(i,ind,:)= w'*L(i,:);
    M (i,ind)  = w;  
    
%     for j = 1:r
%         [junk, ind] = min(dist);
%         dist(ind)   = inf;
%         OL(i,ind,:) = full(w(r)*L(i,:))';
%         M (i,ind)   = w(r); % return this for future calls
%     end
end
fprintf('\n');

if ~interpl
    OL = mean((OL),1); % mean value of a given vertex
else
    fprintf('Averaging proximal vertices (wait...)');
    for i = 1:size(OL,2)
        for j = 1:size(OL,3)
            % average overlapping voxels 
            L(i,j) = sum( OL(:,i,j) ) / length(find(OL(:,i,j))) ;
        end
    end
    fprintf(' ...Done\n');
    OL = L;
end

% normalise and rescale
for i = 1:size(OL,2)
    this = OL(:,i);
    y(:,i)  = S(i,1) + ((S(i,2)-S(i,1))).*(this - min(this))./(max(this) - min(this));
end

y(isnan(y)) = 0;
y  = full(y);
%bigimg;view(270,0);


% spm mesh smoothing
fprintf('Smoothing overlay...\n');
for i = 1:ntime
    y(:,i) = spm_mesh_smooth(mesh, double(y(:,i)), 4);
end

% close image so can reopen with subplots
if num == 2;
    close
    f  = figure;
    set(f, 'Position', [100, 100, 2000, 1000])
    h1 = subplot(121);
    h2 = subplot(122);
else
    bigimg;view(270,0);
    f = gcf;
end

% MAKE THE GRAPH / VIDEO
D        = 'headlight';
cam      = camlight(D);
vidObj   = VideoWriter(fpath,'MPEG-4');
set(vidObj,'Quality',100);
set(vidObj,'FrameRate',size(y,2)/(10));
open(vidObj);

for i = 1:ntime
    
    if i > 1; fprintf(repmat('\b',[1 length(str)])); end
    str = sprintf('building: %d of %d\n',i,ntime);
    fprintf(str);
    camlight(cam); 
    
    switch num
        case 2
            plot(h1,gifti(mesh));
            hh       = get(h1,'children');
            set(hh(end),'FaceVertexCData',y(:,i), 'FaceColor','interp');    
            shading interp
            view(270,0);
            caxis([min(S(:,1)) max(S(:,2))]);


            plot(h2,gifti(mesh));
            hh       = get(h2,'children');
            set(hh(3),'FaceVertexCData',y(:,i), 'FaceColor','interp');    
            shading interp
            view(90,0);
            caxis([min(S(:,1)) max(S(:,2))]);

        
        case 1
            hh = get(gca,'children');
            set(hh(end),'FaceVertexCData',y(:,i), 'FaceColor','interp');
            caxis([min(S(:,1)) max(S(:,2))]);
            shading interp
    end
    
    try
        tt = title(num2str(tv(i)),'fontsize',20);
        P = get(tt,'Position') ;
        P = P/max(P(:));
        set(tt,'Position',[P(1) P(2)+70 P(3)]) ;
    end
    
    set(findall(gca, 'type', 'text'), 'visible', 'on');
    
    if colbar
        colorbar
    end
    drawnow;
            
              

    currFrame = getframe(f);
    writeVideo(vidObj,currFrame);
end
close(vidObj);


    
end
















% Notes / Workings
%---------------------------------------------------
    %rotations - because x is orientated backward?
%     t  = 90;
%     Rx = [ 1       0       0      ;
%            0       cos(t) -sin(t) ;
%            0       sin(t)  cos(t) ];
%     Ry = [ cos(t)  0      sin(t)  ;
%            0       1      0       ;
%           -sin(t)  0      cos(t)  ];
%     Rz = [ cos(t) -sin(t) 0       ;
%            sin(t)  cos(t) 0       ;
%            0       0      1       ];
   %M = (Rx*(M'))';
   %M = (Ry*(M'))';
   %M = (Rz*(M'))';

   
