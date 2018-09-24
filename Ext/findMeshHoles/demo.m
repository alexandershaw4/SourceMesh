%% DEMO
function demo()
% load sample mesh - created from circleMesh.m
load circleMeshSample.mat

xy = p';
faces = t(1:3,:)';

% add 3rd dimension (z)
z = ones(size(xy,1),1);
vertices = [xy,z];

figure; trimesh(faces,vertices(:,1),vertices(:,2),vertices(:,3)); axis equal;
title('Original'), xlabel('x'), ylabel('y'), zlabel('z')

% create hole
holeSize = (max(p(1,:))-min(p(2,1)))/4;
[vertices,faces] = createHole(vertices,faces,holeSize);
figure; trimesh(faces,vertices(:,1),vertices(:,2),vertices(:,3)); axis equal;
title('Add Hole'), xlabel('x'), ylabel('y'), zlabel('z')

% find hole
holeCellArray = findTriMeshHoles(faces,vertices);

% view holes
figure; trimesh(faces,vertices(:,1),vertices(:,2),vertices(:,3));
title('Identify Holes'); xlabel('x');ylabel('y'); hold on; axis equal;
for i = 1:length(holeCellArray)
    hole = holeCellArray{i};
    line(vertices(hole,1),vertices(hole,2),vertices(hole,3),'Color','r')
end
end

function [vertices,faces] = createHole(vertices,faces,holeSize)
% create hole in circle mesh

dimens = size(vertices,2);

vertices_norm = vertices - repmat(mean(vertices),size(vertices,1),1);
temp3 = vertices_norm(faces,:);
matrix3 = reshape(temp3,size(faces,1),dimens,3);
faceCenters = sum(matrix3,2)./3;
dist2 = sqrt(sum(faceCenters.^2,3));
faces = faces(dist2 >= holeSize,:);

end