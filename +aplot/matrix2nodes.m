function [node1,node2,strng] = matrix2nodes(A,pos)
% Write node & edge files for the AAL90 atlas
% Also returns node-to-node coordinates for the matrix specified.
%
% Input is the n-by-n connectivity matrix
% Input 2 is the sourcemodel vertices, n-by-3
%
% AS2017



node1 = []; node2 = []; strng = [];
for i = 1:length(A)
    [ix,iy,iv] = find(A(i,:));
    
    if ~isempty(ix)
        conns = max(length(ix),length(iy));
        for nc = 1:conns
            node1 = [node1; pos(i(1),:)];
            node2 = [node2; pos(iy(nc),:)];
            strng = [strng; iv(nc)];
        end
    end
end

end