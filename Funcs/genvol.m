function V = genvol(vertices,values,dims)

% initialise volume
V = zeros(dims);

% make a grid in each direction, with some whitespace at edges
lb = min(vertices)*1.2;
ub = max(vertices)*1.2;

% using common x,y,z lengths will preserve ratios
edges = [ min(lb) max(ub) ];
x = linspace( edges(1),edges(2),dims(1) );
y = linspace( edges(1),edges(2),dims(2) );
z = linspace( edges(1),edges(2),dims(3) );
v = [x' y' z'];

% from our new grid, compute closest points of vertices of interest
% nb. indices of new grid are the volume indices

%D0 = cdist(vertices,v);

% now loop to find clsest points
for i = 1:3
    for j = 1:length(vertices)
        v0(j,i) = findthenearest( vertices(j,i), v(:,i) );
    end
end

% now embed these positions in the volume
for i = 1:length(vertices)
    V( v0(i,1) , v0(i,2) , v0(i,3) ) = values(i);
end