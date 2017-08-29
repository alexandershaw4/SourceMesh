
S   = load('T-tests');
Sub = 2;
t   = S.t(Sub,:);
p   = S.p(Sub,:);

% For top n-points:
%----------------------------------------------------------------------
% [v,i] = maxpoints(abs(t),10);
% N     = zeros(90,1);
% N(i)  = 1;
% t     = double(t);
% atemplate('overlay',t,'nodes',N,'labels');
% bigimg;


% For significant points only:
%----------------------------------------------------------------------
t = double(mean(S.t));

%i = find(mean(S.p)<.05);
i = find(t>2);
N = zeros(90,1);
N(i) = 1;
atemplate('overlay',t,'nodes',N,'labels');
bigimg;

% For Average over subjects
%----------------------------------------------------------------------
muT = mean(double(S.t),1);
N = zeros(90,1);
N([11 12 15 16 17 18 79 80 81 82 85 86]) = 1;

atemplate('overlay',muT,'nodes',N,'labels');
bigimg;

% For MMN-ERP network [A1, STG, IFG] overlay with network and labels
%----------------------------------------------------------------------
S   = load('T-tests');
t   = mean(S.t);
t   = double(t);
load labels

N     = zeros(90,1);
Gr    = [11 12 81 82 79 80];
N(Gr) = 1;
    
% %    L R L R L R
% %    i i A A S S
% e = [0 0 0 0 2 0; % i L
%      0 0 0 0 0 2; % i R
%      0 0 0 0 1 0; % A L
%      0 0 0 0 0 1; % A R
%      1 0 2 0 0 0; % S L
%      0 1 0 2 0 0];% S R
% l = labels(Gr);
% [nodes,edges] = write_subset_nodes(l,e);

atemplate('overlay',t,'nodes',N,'labels');
