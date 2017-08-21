

[v,i] = maxpoints(abs(t),10);
N     = zeros(90,1);
N(i)  = 1;
t     = double(t);

atemplate('overlay',t,'nodes',N,'labels');

