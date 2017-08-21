
S   = load('T-tests');
Sub = 2;
t   = S.t(Sub,:);
p   = S.p(Sub,:);

% Top n-points:
%----------------
% [v,i] = maxpoints(abs(t),10);
% N     = zeros(90,1);
% N(i)  = 1;
% t     = double(t);

% Signif points:
%----------------
t = double(t);
i = find(p<.05);
N = zeros(90,1);
N(i) = 1;


atemplate('overlay',t,'nodes',N,'labels');


muT = mean(double(S.t),1);
atemplate('overlay',muT);
bigimg;
