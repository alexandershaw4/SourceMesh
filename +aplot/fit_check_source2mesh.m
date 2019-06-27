function pos = fit_check_source2mesh(pos,g)
% check box bounds of sourcemodel and mesh are matched
%
% pos = nx3 source model, g = mesh gifti/struct

% ensure sourcemodel (pos) is around same scale as mesh boundaries
m = min(g.vertices);% *1.1;
M = max(g.vertices);% *1.1;

% refit box boundaries
V        = pos - repmat(spherefit(pos),[size(pos,1),1]);
V(:,1)   = m(1) + ((M(1)-m(1))).*(V(:,1) - min(V(:,1)))./(max(V(:,1)) - min(V(:,1)));
V(:,2)   = m(2) + ((M(2)-m(2))).*(V(:,2) - min(V(:,2)))./(max(V(:,2)) - min(V(:,2)));
V(:,3)   = m(3) + ((M(3)-m(3))).*(V(:,3) - min(V(:,3)))./(max(V(:,3)) - min(V(:,3)));
pos      = V;

end
