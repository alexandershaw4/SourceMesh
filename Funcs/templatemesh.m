function templatemesh(A,labelflag)
% Overlay a 90x90 AAL connectivity matrix on a template brain natively in
% matlab 
%
%
% AS17


% Surface
mesh = read_nv();
meshmesh(mesh); 
hold on;


% AAL Connectivity
if nargin > 0

    % Edges
    %-------------------------------------------------------------
    [node1,node2,strng] = conmat2nodes(A);
    RGB = makecolbar(strng);
    
    % LineWidth (scaled) for strength 
    if any(strng)
        R = [min(strng),max(strng)];
        S = ( strng - R(1) ) + 1e-3;
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
    
    % Add Labels
    if nargin > 1 && labelflag == 1;
        addlabels(A);
    end

    drawnow;
    
    
    % Nodes (of edges only)
    %-------------------------------------------------------------
    for i = 1:size(node1,1)
        scatter3(node1(i,1),node1(i,2),node1(i,3),'filled','k');
        scatter3(node2(i,1),node2(i,2),node2(i,3),'filled','k');
    end    
    
end

end

function RGB = makecolbar(I)
% Register colorbar values to our T-vector

Colors   = jet;
NoColors = length(Colors);

Ireduced = (I-min(I))/(max(I)-min(I))*(NoColors-1)+1;
RGB      = interp1(1:NoColors,Colors,Ireduced);

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
off = 1.1;

% add these to plot with offset
%------------------------------------
for i = 1:length(AN)
    t(i) = text(v(AN(i),1)+off,+v(AN(i),2)+off,v(AN(i),3)+off,labels{AN(i)});
end
set(t,'Fontsize',7)

end
