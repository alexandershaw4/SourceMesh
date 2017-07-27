function templatemesh(A)
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
    S = (strng/max(strng))*7;
    R = [min(strng),max(strng)];
    S = ( strng - R(1) ) + 1e-3;
    
    for i = 1:size(node1,1)
        line([node1(i,1),node2(i,1)],...
             [node1(i,2),node2(i,2)],...
             [node1(i,3),node2(i,3)],...
             'LineWidth',S(i),'Color',[RGB(i,:)]);
    end
    set(gcf,'DefaultAxesColorOrder',RGB)
    colorbar
    caxis([min(strng),max(strng)])

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
% V = vsmooth(g.vertices, g.faces, .06);
% g.vertices = V;

h = plot(gifti(g));
C = [.5 .5 .5];

set(h,'FaceColor',[C]); box off;
grid off;  set(h,'EdgeColor','none');
alpha(.5); set(gca,'visible','off');

end
