function templatemesh(A);

% Surface
mesh = read_nv();
meshmesh(mesh); 
hold on;

if nargin > 0
    % Nodes
    load AAL_SOURCEMOD
    N = template_sourcemodel.pos;
    for i = 1:90
        scatter3(N(i,1),N(i,2),N(i,3),'filled','r');
    end
    drawnow;
    
    % Edges
    [node1,node2,strng] = conmat2nodes(A);
    for i = 1:size(node1,1)
        line([node1(i,1),node2(i,1)],[node1(i,2),node2(i,2)],[node1(i,3),node2(i,3)]);
    end
    drawnow;
end
end



function meshmesh(g)


% f = g.faces;
% v = g.vertices;
% h = patch('faces',f,'vertices',v);
h = plot(gifti(g));
C = .5;

set(h,'FaceColor',[C C C]);
box off;
grid off; 
set(h,'EdgeColor','none')
alpha(.3);
set(gca,'visible','off');

end
