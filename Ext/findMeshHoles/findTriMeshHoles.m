function [holeCellArray,bounding_triangles,holeLengths] = findTriMeshHoles(faces,vertices)
% Finds holes in a triangular mesh 
% Note: Does not work if a hole shares more than one vertex with other holes 
% Input: 
% faces = M x 3 
% vertices = N x 3 (optional if you want the hole lengths)

% Output: 
% holeCellArray = P x 1 cell array containing a list of holes, which are 
% traced in consecutive order (list of scalar indices) 
% bounding_triangles = Q x 3 list of faces that contain a bounding edge (does 
% not contain triangles that only has a single bounding vertex) 
% holeLengths = P x 1 vector containing the perimeter of each hole

if nargout > 2 && nargin < 2
    error('Please input vertices if you want the hole lengths')
end

% sort each row
faces = sort(faces,2);

% face error checking (ensure no faces contain an index more than once)
sel = faces(:,1) == faces(:,2) | faces(:,2) == faces(:,3) | faces(:,1) == faces(:,3);
if sum(sel) > 0
    warning('Triangle refers to a vertex more than once... Consider fixing your data...')
    faces = faces(~sel,:);
end

% get edges from faces
edge1 = faces(:,1:2);
edge2 = faces(:,[1,3]);
edge3 = faces(:,2:3);
edges = vertcat(edge1,edge2,edge3);

% find edges that have no replicates (find boundary edges)
edges = sortrows(edges);
vec = diff(edges); % diff btw rows
vec = any(vec,2); % unique rows
vec1 = [1;vec]; % doesn't match with previous row
vec2 = [vec;1]; % doesn't match with next row
vec = vec1&vec2; % doesn't match previous and next row
bounding_edges = edges(vec,:); % unique edges

% bounding triangles - triangles that contain a bounding edge 
% Note: you can modify to get triangles that contain a bounding vertex
match1 = ismember(edge1,bounding_edges,'rows');
match2 = ismember(edge2,bounding_edges,'rows');
match3 = ismember(edge3,bounding_edges,'rows');
bounding_triangles = faces(match1|match2|match3,:);

if nargin > 1 && nargout > 1
    % length of each edge
    edgeVertex1 = vertices(bounding_edges(:,1),:);
    edgeVertex2 = vertices(bounding_edges(:,2),:);
    edge_lengths = sqrt(sum((edgeVertex1-edgeVertex2).^2,2));
    [holeCellArray,holeLengths]=getTraces(bounding_edges,edge_lengths);
else
    holeCellArray=getTraces(bounding_edges);
end


end

function [tracesDB,lengthDB,completeDB]=getTraces(bounding_edges,edge_lengths)
% Connects the unordered edges to form holes that list the edges in
% consecutive order.
% The number of formed holes from the edges determines the size of tracesDB

% INPUT
% bounding_edges = N x 2 list of connected border vertices
% edge_lengths   = N x 1 the corresponding length of the connected vertices

% OUTPUT
%   tracesDB = cell array containing a list of holes, which are traced in
%               consecutive order
%   lengthDB = 1-D vector containing the perimeter of each hole that
%                   corresponds to tracesDB
%   completeDB = 1-D vector that indicates if corresponding hole in tracesDB is completely traced
%                 (technically should be all 1's) - for testing purposes

if nargin < 2
    edge_lengths = zeros(size(bounding_edges,1),1);
end

tracesDB = {};
lengthDB = [];
completeDB = [];
trace12 = [];
idx = 1;
while size(bounding_edges,1) > 0
    trace = bounding_edges(1,:)'; % add first trace
    bounding_edges = bounding_edges(2:end,:);
    total_length = edge_lengths(1); % add first trace length 
    edge_lengths = edge_lengths(2:end);
    flip_flag = 0;
    complete_flag = 1;
    while trace(1) ~= trace(end)
        [row_idx,col_idx] = find(bounding_edges == trace(end));
        col_idx = -col_idx+3;
        % if more than one bounding edge connects to the end of the
        % trace then trace in the other direction (occurs when only one
        % triangle separates two holes so both holes share a single vertex)
        if length(row_idx) > 1 || length(row_idx) < 1 
            if flip_flag == 0
                trace = flip(trace);
                flip_flag = 1;
                continue;
            else
                complete_flag = 0;
    %             disp('Incomplete hole tracing... May need to make more conditions')
                break;
            end
        end
        trace = [trace;bounding_edges(row_idx,col_idx)];
        total_length = total_length+edge_lengths(row_idx);
        bounding_edges(row_idx,:) = [];
        edge_lengths(row_idx) = [];
    end
    tracesDB = [tracesDB;{trace}]; 
    lengthDB = [lengthDB; total_length]; 
    completeDB = [completeDB;complete_flag]; 
    trace12 = [trace12;trace(1),trace(end)]; 
    idx = idx+1;
end
if sum(~completeDB) > 0
    warning('Not all holes are completely traced due to holes sharing more than one vertex with other holes')
end
end
