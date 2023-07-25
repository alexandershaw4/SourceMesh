function quickscatter3(v,varargin)

if nargin == 1

    scatter3(v(:,1),v(:,2),v(:,3),100);

else
    scatter3(v(:,1),v(:,2),v(:,3),100);hold on;
    for i = 1:length(varargin)
        vx = varargin{i};
        scatter3(vx(:,1),vx(:,2),vx(:,3),100);hold on;
    end

end