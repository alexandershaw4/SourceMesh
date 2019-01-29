function newpos = fixmesh(g,pos)
% plot as transparent grey gifti surface
%
% AS

v = g.vertices;
v = v - repmat(spherefit(v),[size(v,1),1]); % Centre on ~0
g.vertices=v;

% Centre on ~0
pos = pos - repmat(spherefit(pos),[size(pos,1),1]);

for i = 1:length(pos)
    this  = pos(i,:);
    [t,I] = maxpoints(cdist(v,this),1,'max');
    newpos(i,:) = v(I,:);
end

end