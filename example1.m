% Synthetic toy point cloud
N = 10000;
theta = 2*pi*rand(N,1);
phi = acos(2*rand(N,1)-1);
r = 1 + 0.05*randn(N,1);  % noisy radius

X = [r.*sin(phi).*cos(theta), ...
     r.*sin(phi).*sin(theta), ...
     r.*cos(phi)];

% Colour by spatial position (so nearby points have similar colours)
C = 0.5 + 0.5*[X(:,1), X(:,2), X(:,3)];
C = max(min(C,1),0);  % clamp to [0,1]

% Fit variational Bayesian Gaussian splats
K = 64; iters = 20;
[M, hist] = vbgs_fit_batch(X, C, K, iters);

% Predict colour for some query points
[c_hat, r] = vbgs_query_color(M, X(1:10,:));
disp(c_hat(1:5,:))

% Visualise components (means)
scatter3(X(:,1), X(:,2), X(:,3), 5, C, 'filled'); hold on;
mu = cat(2,M.spat_post.mu)';
scatter3(mu(:,1), mu(:,2), mu(:,3), 100, 'k', 'filled');
axis equal; title('VBGS synthetic sphere');