function data = connections(data,A,colbar,write,fname,netcmap)
% Network (Node & Edges) plotter function
%
%
pos = data.sourcemodel.pos;

% Read edge/node files if string
%--------------------------------------------------------------------------
if ischar(A)
    [fp,fn,fe]  = fileparts(A);
    [edge,node] = aplot.rw_edgenode(fn);
    A           = edge;
    
    if isfield(data,'template')
        if isfield(data.template,'model')
           fprintf('Doing atlas registration\n');
           i.template = 1;
           i.model    = data.template.model;
           i.labels   = data.template.labels;
           i.A        = A;
           data.sourcemodel.pos = node(:,1:3);
           [data,i]   = aplot.sort_template(data,i);
           A          = i.A;
        end
    end
end

A(isnan(A)) = 0;
A(isinf(A)) = 0;

% rescale network positions inside box boundaries of mesh
% (i thought meshmesh had already done this?)
if ~isempty(data.mesh.h)
    bounds = [min(data.mesh.vertices); max(data.mesh.vertices)];
    offset = 0.99;
    for ip = 1:3
        pos(:,ip) = bounds(1,ip) + ((bounds(2,ip)-bounds(1,ip))) .* ...
                    (pos(:,ip) - min(pos(:,ip)))./(max(pos(:,ip)) - min(pos(:,ip)));
        pos(:,ip) = pos(:,ip)*offset;
    end

    % redirect to closest mesh point (vertex)
    for ip = 1:length(pos)
        [~,this]  = min(cdist(pos(ip,:),data.mesh.vertices));
        pos(ip,:) = data.mesh.vertices(this,:);
    end
end

% Edges
%--------------------------------------------------------------------------
[node1,node2,strng] = aplot.matrix2nodes(A,pos);

% place both signed absmax value in overlay so that colorbar is symmetrical
strng2 = [strng; -max(abs(strng)); max(abs(strng))];
RGB    = aplot.makecolbar(strng2,netcmap);

% LineWidth (scaled) for strength
if any(strng)
    R = [-max(abs(strng)),max(abs(strng))];
    S = ( abs(strng) - R(1) ) + 1e-3;
    
    % If all edges same value, make thicker
    if  max(S(:)) == 1e-3; 
        S = 3*ones(size(S)); 
    end
else
    S = [0 0];
end

% If too few strengths, just use red edges
%--------------------------------------------------------------------------
LimC = 1;
if all(all(isnan(RGB)))
    RGB  = repmat([1 0 0],[size(RGB,1) 1]);
    LimC = 0;
end

data.network.edge = A;
data.network.node = pos;
data.network.RGB  = RGB;
data.network.tofrom.node1 = node1;
data.network.tofrom.node2 = node2;

% force girth boundaries if stable
if ~any(isnan( (S - min(S)) ./ (max(S) - min(S)) ))
    S = 0.1 + (3 - 0) .* (S - min(S)) ./ (max(S) - min(S));
end

% Paint edges
%--------------------------------------------------------------------------
for i = 1:size(node1,1)
    l0(i)=line([node1(i,1),node2(i,1)],...
        [node1(i,2),node2(i,2)],...
        [node1(i,3),node2(i,3)],...
        'LineWidth',S(i),'Color',[RGB(i,:)]);
end

% Set colorbar only if there are valid edges
%--------------------------------------------------------------------------
if any(i) && colbar
    set(gcf,'DefaultAxesColorOrder',RGB)
    set(gcf,'Colormap',RGB)
    if colbar
        %colormap(jet)
        %colorbar
        drawnow; pause(.5);
        a1  = gca;
        axb = axes('position', get(a1, 'position'));
        set(axb,'visible','off')
        axes(axb);
        %set(a1,'DefaultAxesColorOrder',RGB)
        set(gcf,'Colormap',RGB)
        
        if any(any(netcmap ~= 0)); 
                    colormap(netcmap);
        else;       colormap(jet);
        end
        
        colorbar('peer',a1,'South');
    end
end
if LimC && colbar
    axes(a1);
    caxis(R);
end

drawnow;


% Nodes (of edges only)
%--------------------------------------------------------------------------
% hold on;
% for i = 1:size(node1,1)
%     scatter3(node1(i,1),node1(i,2),node1(i,3),'filled','k');
%     scatter3(node2(i,1),node2(i,2),node2(i,3),'filled','k');
% end

drawnow;

if write;
   fprintf('Writing network: .edge & .node files\n');
   conmat2nodes(A,fname,'sourcemodel',pos);
end


end