function data = drawnodes(data, N)
% Node plotter. N = (90,1) with 1s for nodes to plot and 0s to ignore.
%
% 
hold on;
pos = data.sourcemodel.pos;
%v   = pos*0.9;

if isfield(data.sourcemodel,'vi')
    % if pos is cell, it's because we've passed both a vertex set and a
    % vector that describes which vertex belongs to which roi
    v  = pos;
    vi = data.sourcemodel.vi;
    n  = unique(vi);
    for i = 1:length(n)
        these = find(vi==n(i));
        roi(i,:) = spherefit(v(these,:));
    end
    pos = roi;
    for ip = 1:length(pos)
        [~,this]  = min(cdist(pos(ip,:),data.mesh.vertices));
        pos(ip,:) = data.mesh.vertices(this,:);
    end
end


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
    
elseif iscell(N)
    % pass a cell of length sm with color strings
    % e.g. N=cell(90,1);
    %      N{1} = 'r';
    ForPlot = [];
    for j = 1:length(N)
        if any(N{j})
            v0 = v(j,:);
            ForPlot = [ForPlot v0];
            scatter3(v0(1),v0(2),v0(3),150,N{j},'filled');
        end
    end
    
else
    ForPlot = v(find(N),:);
    %s       = find(N);
    s = ones(length(find(N)),1)*150;
    for i   = 1:size(ForPlot,1)
        col = 'r';
        scatter3(ForPlot(i,1),ForPlot(i,2),ForPlot(i,3),s(i),'r','filled');
    end
end
%RGB = makecolbar(ForPlot);
%set(gcf,'DefaultAxesColorOrder',RGB); jet;
colorbar

data.drawnodes.data = ForPlot;

end