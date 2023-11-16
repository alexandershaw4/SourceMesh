function v0 = alignclouds(v0,v)
% 
% v0 = alignclouds(v0,v)
%
%

v0 = centralise(v0) ;
v0 = boxbound(v,v0);

M = ones(10,1);

f = @(M) p2pmap(v,v0*reshape(M(1:9),[3 3])*M(10));

%[X, fX, i] = PR_minimize(M, f, 1)

end

function e = p2pmap(v,v0)

[IDX, D] = knnsearch(v,v0);
e = norm(D);

end

function g = grid3(b,n)
for i = 1:3    
    g(:,i) = linspace(b(1,i),b(2,i),n);
end
end

function v = centralise(v)

v  = v - repmat(spherefit(v),[size(v,1),1]);

end

function [v0,b] = boxbound(v,v0)

b  = [min(v);max(v)];

for i = 1:3
    v0(:,i) = rescale(v0(:,i),b(1,i),b(2,i));
end

end

function [vin,Ed] = dosphere(v)

    Bnd = [min(min(v)); max(max(v))]*1.1;
    Ed  = cdist(v,spherefit(v));
    v   = v ./ [Ed/3 Ed/3 Ed/3];
    for i = 1:3
        vin(:,i) = rescale(v(:,i),Bnd(1),Bnd(2));
    end

end

function [Centre,A,B] = spherefit(X)
% Fit sphere to centre of vertices, return centre points
%
%

A =  [mean(X(:,1).*(X(:,1)-mean(X(:,1)))), ...
    2*mean(X(:,1).*(X(:,2)-mean(X(:,2)))), ...
    2*mean(X(:,1).*(X(:,3)-mean(X(:,3)))); ...
    0, ...
    mean(X(:,2).*(X(:,2)-mean(X(:,2)))), ...
    2*mean(X(:,2).*(X(:,3)-mean(X(:,3)))); ...
    0, ...
    0, ...
    mean(X(:,3).*(X(:,3)-mean(X(:,3))))];
A = A+A.';
B = [mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,1)-mean(X(:,1))));...
     mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,2)-mean(X(:,2))));...
     mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,3)-mean(X(:,3))))];
Centre=(A\B).';
end

function [e0,g] = dfdx(f,x)
% simple 1-step finite difference routine for compute partial gradients

e0 = f(x);
k  = exp(-8);

for i  = 1:length(x)
    dx    = x;
    dx(i) = dx(i) + k;
    g(i)  = (f(dx) - e0) ./ k;
end
g=g';

end