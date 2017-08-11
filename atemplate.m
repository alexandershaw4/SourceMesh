function atemplate(varargin)
% add AAL networks and overlays to a smoothed brain
%
% paired inputs:
%
% 'overlay', L (1x90)
% 'network', A (90x90)
% 'tracks' , tracks, header (as read with trk_read^ )
% 'labels'
%
% usages:
%
%   atemplate('labels');     template mesh with AAL labels
%   atemplate('overlay',L);  template mesh with overlay
%   atemplate('network',A);  template mesh nodes & edges
%   atemplate('overlay',L,'network',A,'labels'); with overlay,network and labels
%   atemplate('tracks',tracks,header); plot tracks loaded with trk_read
%
%
% ^trk_read requires along-tract-stats toolbox
%
% AS17

pmesh  = 1;
labels = 0;
for i  = 1:length(varargin)
    if strcmp(varargin{i},'overlay'); L = varargin{i+1}; end
    if strcmp(varargin{i},'network'); A = varargin{i+1}; end
    if strcmp(varargin{i},'tracks');  T = varargin{i+1}; H = varargin{i+2}; end
    if strcmp(varargin{i},'labels');  labels = 1;  end
    if strcmp(varargin{i},'nosurf');  pmesh  = 0;  end
    if strcmp(varargin{i},'nodes');   N = varargin{i+1}; end
end


% Plot Surface
mesh = read_nv();
hold on;

if     pmesh && ~exist('T','var');
       meshmesh(mesh);
elseif pmesh
       meshmesh(mesh,.2);
end

try A; connections(A);        end % draw edges and edge-connected nodes
try L; overlay(mesh,L);       end % find closest vertices and overlay
try T; drawtracks(T,H,mesh);  end % draw dti tracks loaded with trk_read
try N; drawnodes(N);          end

if labels; 
    if     exist('A','var'); addlabels(A);
    elseif exist('N','var'); addlabels(diag(N));
    else   addlabels(ones(90,90));
    end
end


end



function connections(A)

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
if any(i)
    set(gcf,'DefaultAxesColorOrder',RGB)
    colorbar
end
if LimC;
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

hold on;
load([fileparts(which('conmat2nodes')),'/AAL_SOURCEMOD.mat']);
v = template_sourcemodel.pos;

ForPlot = v(find(N),:);

for i = 1:length(ForPlot)
    scatter3(ForPlot(i,1),ForPlot(i,2),ForPlot(i,3),'filled','r');
end


end

function RGB = makecolbar(I)
% Register colorbar values to our T-vector

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

% this time draw the tracks
for iTrk = 1:length(tracks)
    matrix = tracks(iTrk).matrix;
    matrix(any(isnan(matrix(:,1:3)),2),:) = [];
    
    matrix(:,2) = matrix(:,2)*-1;
    M           = matrix - repmat(Centre,[size(matrix,1),1]); % centre
    M(:,2)      = M(:,2) + pullback;                          % pullback
    M(:,3)      = M(:,3) + pullup;                            % pull up

    h = patch([M(:,1)' nan], [M(:,2)' nan], [M(:,3)' nan], 0);
    cdata = [(0:(size(matrix, 1)-1))/(maxpts) nan];
    set(h,'cdata', cdata, 'edgecolor','interp','facecolor','none');
end

h = get(gcf,'Children');
set(h,'visible','off');

end


function overlay(mesh,L);

% overlay option:
interpl = 0; % interp shading between nodes or just use max value

% Overlay
load('AAL_SOURCEMOD');
v  = template_sourcemodel.pos;
x  = v(:,1);
mv = mesh.vertices;
nv = length(mv);
OL = sparse(90,nv);
r  = 1500;          % radius func
w  = linspace(.1,1,r);
w  = fliplr(w); 
M  = zeros( length(x), nv);

% if is same verts as mri, just rescale & overlay
if length(L) == nv && nargin == 2
    fprintf('Overlay size matches mri!\n');
    hh = get(gca,'children');
    L  = L(:);
    y  = S(1) + ((S(2)-S(1))).*(L - min(L))./(max(L) - min(L));
    y  = y(:);

    set(hh(end),'FaceVertexCData',y, 'FaceColor','interp');
    shading interp;
    
else
    S  = [min(L(:)),max(L(:))];
    % otherwise find closest points (assume both in mm)
    fprintf('Determining closest points between AAL & template vertices (overlay)\n');
    for i = 1:length(x)

        % reporting
        if i > 1; fprintf(repmat('\b',[size(str)])); end
        str = sprintf('%d/%d',i,(length(x)));
        fprintf(str);    

        % find closest point in cortical mesh
        dist  = cdist(mv,v(i,:));
        for j = 1:r
            [junk, ind] = min(dist);
            dist(ind)   = inf;
            OL(i,ind)   = w(r)*L(i);
            M (i,ind)   = w(r); % return this for future calls
        end
    end
    
    if ~interpl
        OL = max(OL)'; % max value of a given vertex
    else
        for i = 1:size(OL,2)
           % average overlapping voxels 
           L(i) = sum( OL(:,i) ) / length(find(OL(:,i))) ;
        end
        OL = L;
        OL = full(L);
    end
    
    % normalise and rescale
    y  = S(1) + ((S(2)-S(1))).*(OL - min(OL))./(max(OL) - min(OL));
    y  = y(:);
    hh = get(gca,'children');
    
    set(hh(end),'FaceVertexCData',y, 'FaceColor','interp');
    shading interp
    colorbar
end

end




function meshmesh(g,a)
% plot as transparent grey gifti surface
%
% AS

% Smooth brain?
% V = vsmooth(g.vertices, g.faces, .03);
% g.vertices = V;

if nargin == 1;
    a = .5;
end

h = plot(gifti(g));
C = [.5 .5 .5];

set(h,'FaceColor',[C]); box off;
grid off;  set(h,'EdgeColor','none');
alpha(a); set(gca,'visible','off');

h = get(gcf,'Children');
set(h(end),'visible','off');


end

function addlabels(V)
% find indices of nodes->labels->plot
%
%

load('labels');
labels = strrep(labels,'_',' ');

load('AAL_SOURCEMOD');

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
v   = template_sourcemodel.pos;
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
set(t,'Fontsize',7)

end

function Centre = spherefit(X)

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

   
