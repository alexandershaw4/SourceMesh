function template(varargin)
% add AAL networks and overlays to a smoothed brain
%
% paired inputs:
%
% 'overlay', L (1x90)
% 'network', A (90x90)
% 'labels'
%
% usages:
%
%   template('labels');     template mesh with AAL labels
%   template('overlay',L);  template mesh with overlay
%   template('network',A);  template mesh nodes & edges
%   template('overlay',L,'network',A,'labels'); with overlay,network and labels
%
% AS17

for i = 1:length(varargin)
    if strcmp(varargin{i},'overlay'); L = varargin{i+1}; end
    if strcmp(varargin{i},'network'); A = varargin{i+1}; end
    if strcmp(varargin{i},'tracks');  T = varargin{i+1}; H = varargin{i+2}; end
    if strcmp(varargin{i},'labels');  labels = 1; else labels = 0; end
end

% Plot Surface
mesh = read_nv();
meshmesh(mesh); 
hold on;

try A; connections(A);   end
try L; overlay(mesh,L);  end
try T; drawtracks(T,H,mesh);  end

if labels; 
    try   A; addlabels(A); 
    catch    addlabels(ones(90,90));
    end
end


end



function connections(A)
% AAL Connectivity

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
for i = 1:size(node1,1)
    scatter3(node1(i,1),node1(i,2),node1(i,3),'filled','k');
    scatter3(node2(i,1),node2(i,2),node2(i,3),'filled','k');
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
hold on; clc;

A   = header.vox_to_ras;
A   = A(1:3,4)';
All = [];

for iTrk = 1:length(tracks)
    if iTrk > 1; fprintf(repmat('\b',size(str))); end
    str = sprintf('Compiling %d of %d\n',iTrk,length(tracks));
    fprintf(str);
    
    matrix = tracks(iTrk).matrix;
    matrix(any(isnan(matrix(:,1:3)),2),:) = [];
    matrix = matrix./(repmat(A,[size(matrix,1),1]));

    All = [All ; matrix];
    %line(matrix(:,1), matrix(:,2), matrix(:,3));
end

% centre on 0
mx = median(All(:,1));
my = median(All(:,2));
mz = median(All(:,3));
dAll = All + repmat([mx my mz],[size(All,1),1]);

% expand to template
M(1,:) = min(mesh.vertices);
M(2,:) = max(mesh.vertices);
T      = dAll;

x = M(1,1) + ((M(2,1)-M(1,1))).*(T(:,1) - min(T(:,1)))./(max(T(:,1)) - min(T(:,1)));
y = M(1,2) + ((M(2,2)-M(1,2))).*(T(:,2) - min(T(:,2)))./(max(T(:,2)) - min(T(:,2)));
z = M(1,3) + ((M(2,3)-M(1,3))).*(T(:,3) - min(T(:,3)))./(max(T(:,3)) - min(T(:,3)));

line(x,y,z);

end


function overlay(mesh,L);
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
        dist  = sum((mv - repmat(v(i, :), size(mv, 1), 1)).^2, 2);
        for j = 1:r
            [junk, ind] = min(dist);
            dist(ind)   = inf;
            OL(i,ind)   = w(r)*L(i);
            M (i,ind)   = w(r); % return this for future calls
        end
    end
    
    % normalise and rescale
    OL = max(OL)';
    y  = S(1) + ((S(2)-S(1))).*(OL - min(OL))./(max(OL) - min(OL));
    y  = y(:);
    hh = get(gca,'children');
    
    set(hh(end),'FaceVertexCData',y, 'FaceColor','interp');
    shading interp
    colorbar
end

end




function meshmesh(g)
% plot as transparent grey gifti surface
%
% AS

% Smooth brain?
% V = vsmooth(g.vertices, g.faces, .03);
% g.vertices = V;

h = plot(gifti(g));
C = [.5 .5 .5];

set(h,'FaceColor',[C]); box off;
grid off;  set(h,'EdgeColor','none');
alpha(.5); set(gca,'visible','off');

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
