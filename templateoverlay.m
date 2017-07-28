function templateoverlay(L)
% Plot smoothed tamplate brain with overlay vector, L.
% L is a vector of length 90
%
%


% Surface
mesh = read_nv();
meshmesh(mesh); 
hold on;

% Overlay
load('AAL_SOURCEMOD');
v  = template_sourcemodel.pos;
x  = v(:,1);
y  = v(:,2);
z  = v(:,3);
mv = mesh.vertices;
nv = length(mv);
OL = sparse(90,nv);
r  = 3000;          % radius func

fprintf('Determining closest points between AAL & template vertices\n');
for i = 1:length(x)
    
    % find closest point in cortical mesh
    dist  = sum((mv - repmat(v(i, :), size(mv, 1), 1)).^2, 2);
    for j = 1:r
        [junk, ind] = min(dist);
        dist(ind)   = inf;
        coord       = mv(ind, :);
        newv(i,:,:) = coord;
        OL(i,ind)   = L(i);
    end
end

OL = mean(OL,1)';
hh = get(gca,'children');
set(hh(end),'FaceVertexCData',OL, 'FaceColor','interp');
shading interp

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