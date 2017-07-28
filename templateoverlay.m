function M = templateoverlay(L)
% Plot smoothed tamplate brain with overlay vector, L.
% L is a vector of length 90
%
% L is mapped to the size of the template by finding the closest points and
% linearly interpreting 
%
% Returns matrix M so that it needn't be recomputed:
% call 1  : M = templateoverlay(L) % computes M and plots
% call 2+ : templateoverlay(L'*M); % much quicker
%
% AS

% Surface
mesh = read_nv();
meshmesh(mesh); 
hold on;

% Overlay
load('AAL_SOURCEMOD');
v  = template_sourcemodel.pos;
x  = v(:,1);
mv = mesh.vertices;
nv = length(mv);
OL = sparse(90,nv);
r  = 3000;          % radius func
w  = linspace(.5,1,r);
M  = zeros( length(x), nv);
S  = [min(L(:)),max(L(:))];

if length(L) == nv
    fprintf('Overlay size matches!\n');
    hh = get(gca,'children');
    L = L(:);
    set(hh(end),'FaceVertexCData',L, 'FaceColor','interp');
    shading interp;
    
else
    fprintf('Determining closest points between AAL & template vertices\n');
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
    OL = mean(OL,1)';
    OL = (OL/min(OL))/(max(OL)-min(OL)) * S;
    hh = get(gca,'children');
    set(hh(end),'FaceVertexCData',OL, 'FaceColor','interp');
    shading interp
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
alpha(.8); set(gca,'visible','off');

end