function [dcloud1] = align_clouds_3d(cloud0,cloud1)

global points0 points0_sub points1 doplot

doplot = 1;

if length(cloud0) < length(cloud1)
    % cloud 0 has to be the denser/bigger cloud
    dcloud = cloud1;
    cloud1 = cloud0;
    cloud0 = dcloud;
end

points0 = cloud0;
points1 = cloud1;

% build rotation model (parameter set)
%-------------------------------------------------------------------------
Rx = eye(3);
Rx(4,:) = 1;
Rx(:,4) = 1;

% approximate same box boundaries
%-------------------------------------------------------------------------
for i = 1:3
    points1(:,i) = min(points0(:,i)) + (max(points0(:,i)) - min(points0(:,i))) ...
        .* ( points1(:,i) - min(points1(:,i)) ) / (max(points1(:,i))-min(points1(:,i)));  
end



% dimensionality reduction
%-------------------------------------------------------------------------
% design a reduced-cloud (subset) matching num polpoints
D = cdist(points0,points1);

for i = 1:size(D,2) 
    % find appropriate mri point for each shape point
    [~,ind] = min(D(:,i));
    closests(i) = ind(1);
end
points0_sub = points0(closests,:);

% initial fit error
%-------------------------------------------------------------------------
e0 = sum( sum(points0_sub - points1) ).^2;
fprintf('Initial (squared) fit error = %d\n',e0);

if doplot == 1
    figure('position',[1073         288         800         673]);
end

% OPTIMISATION
%-------------------------------------------------------------------------
[X, F, i] = PR_minimize(Rx(:), @fitter, 128);
%[X,F]  = fminsearch(@fitter,Rx(:) );
fprintf('Posterior (squared) fit error = %d\n',F(end));

% compute posterior positions
%---------------------------------------------------------
Rx       = reshape(X,[4 4]);  % rotation (affine-like)
pnts     = points1*Rx(1:3,1:3);
%scale0   = Rx(1:3,4);     % scaling
%scale1   = Rx(4,1:3);
%bounds   = [min(points1);max(points1)];

%for i = 1:3
%    LB = bounds(1,i)*scale0(i);
%    UB = bounds(2,i)*scale1(i);
%    pnts(:,i) = LB + (UB - LB) .* ( pnts(:,i) - min(pnts(:,i)) ) / (max(pnts(:,i))-min(pnts(:,i)));  
%end

dcloud1 = pnts;


end


function [e,J] = fitter(Rx)

global points0 points0_sub points1 doplot


% the objective function to minimise - i.e.
%
% argmin: e = sum( mri_points - rotation(shapepoints) ).^2
%
%

% compute location given roation and scale
%--------------------------------------------------------------------
Rx = reshape(Rx,[4 4]);
M  = points1*Rx(1:3,1:3);
% scale0 = Rx(1:3,4);
% scale1 = Rx(4,1:3);
% bounds = [min(points1);max(points1)];
% 
% for i = 1:3
%     LB = bounds(1,i)*scale0(i);
%     UB = bounds(2,i)*scale1(i);
%     M(:,i) = LB + (UB - LB) .* ( M(:,i) - min(M(:,i)) ) / (max(M(:,i))-min(M(:,i)));  
% end

% error to minimise
e = sum( sum(points0_sub - M) ).^2;
       
% display
if doplot == 1
    scatter3(points0(:,1),points0(:,2),points0(:,3),20,'b','filled');hold on;
    scatter3(M(:,1),M(:,2),M(:,3),50,'r','filled'); hold off;
    drawnow;
end

if nargout > 1
    %fprintf('computing gradients\n');
    Rx = Rx(:);
    % compute approx jacobian
    delta = .08;
    for i = 1:length(Rx)
        dRx0    = Rx;
        dRx1    = Rx;
        dRx0(i) = dRx0(i) + delta;
        dRx1(i) = dRx1(i) - delta;

        k0      = fitter(dRx0);
        k1      = fitter(dRx1);
        J(i,:)  = (k0 - k1)/(2*delta);
        
    end
end

end