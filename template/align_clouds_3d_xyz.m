function [dcloud1] = align_clouds_3d_xyz(cloud0,cloud1)

global points0 points0_sub points1 doplot

% same as align_clouds_3d.m but uses parametised rotation matrices to fit
% pitch, roll & yaw, along with 3 scale factors
%
%
% AS

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

%     x y z Mx My Mz mx my mz
Rx = [0 0 0 1  1  1  1  1  1];

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
x = X(1);
y = X(2);
z = X(3);
M = points1;

    Rxx= [ 1       0       0      ;
           0       cos(x) -sin(x) ;
           0       sin(x)  cos(x) ];
    Ryy= [ cos(y)  0      sin(y)  ;
           0       1      0       ;
          -sin(y)  0      cos(y)  ];
    Rzz= [ cos(z) -sin(z) 0       ;
           sin(z)  cos(z) 0       ;
           0       0      1       ];
   
M = (Rxx*(M'))';
M = (Ryy*(M'))';
M = (Rzz*(M'))';

scale0 = Rx(4:9);
scale1 = Rx(4:9);
bounds = [min(points1);max(points1)];

for i = 1:3
    LB = bounds(1,i)*scale0(i);
    UB = bounds(2,i)*scale1(i);
    pnts(:,i) = LB + (UB - LB) .* ( M(:,i) - min(M(:,i)) ) / (max(M(:,i))-min(M(:,i)));
end


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
x = Rx(1);
y = Rx(2);
z = Rx(3);
M = points1;

    Rxx= [ 1       0       0      ;
           0       cos(x) -sin(x) ;
           0       sin(x)  cos(x) ];
    Ryy= [ cos(y)  0      sin(y)  ;
           0       1      0       ;
          -sin(y)  0      cos(y)  ];
    Rzz= [ cos(z) -sin(z) 0       ;
           sin(z)  cos(z) 0       ;
           0       0      1       ];
   
M = (Rxx*(M'))';
M = (Ryy*(M'))';
M = (Rzz*(M'))';

scale0 = Rx(4:9);
scale1 = Rx(4:9);
bounds = [min(points1);max(points1)];

for i = 1:3
    LB = bounds(1,i)*scale0(i);
    UB = bounds(2,i)*scale1(i);
    M(:,i) = LB + (UB - LB) .* ( M(:,i) - min(M(:,i)) ) / (max(M(:,i))-min(M(:,i)));
end

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
    delta = exp(-8);%.08;
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