function data = drawnodes(data, N)
% Node plotter. N = (90,1) with 1s for nodes to plot and 0s to ignore.
%
% 
hold on;
pos = data.sourcemodel.pos;
%v   = pos*0.9;

bounds = [min(data.mesh.vertices); max(data.mesh.vertices)];
offset = 0.99;
for ip = 1:3
    pos(:,ip) = bounds(1,ip) + ((bounds(2,ip)-bounds(1,ip))) .* ...
                (pos(:,ip) - min(pos(:,ip)))./(max(pos(:,ip)) - min(pos(:,ip)));
    pos(:,ip) = pos(:,ip)*offset;
end

% redirect to clseast mesh point (vertex?)
for ip = 1:length(pos)
    [~,this]  = min(cdist(pos(ip,:),data.mesh.vertices));
    pos(ip,:) = data.mesh.vertices(this,:);
end

v = pos;

if size(N,1) > 1 && size(N,2) > 1
    cols = {'r' 'm','y','g','c','b'};
    if size(size(N,2)) == 90
        N = N';
    end
    
    for j = 1:size(N,2)
        ForPlot = v(find(N(:,j)),:) + (1e-2 * (2*j) ) ;
        s       = find(N);
        col     = cols{j};
        for i   = 1:length(ForPlot)
            scatter3(ForPlot(i,1),ForPlot(i,2),ForPlot(i,3),70,col,'filled',...
                'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6);        hold on;
        end
    end
    
else
    ForPlot = v(find(N),:);
    %s       = find(N);
    s = ones(length(find(N)),1)*40;
    for i   = 1:length(ForPlot)
        col = 'r';
        scatter3(ForPlot(i,1),ForPlot(i,2),ForPlot(i,3),s(i),'r','filled');
    end
end
%RGB = makecolbar(ForPlot);
%set(gcf,'DefaultAxesColorOrder',RGB); jet;
colorbar

data.drawnodes.data = ForPlot;

end