function [C,verts] = get_roi_centres(pos,all_roi_tissueindex)
% Find centre points of rois
%
%
v   = pos;
roi = all_roi_tissueindex;

i   = unique(roi);
i(find(i==0))=[];

fprintf('Finding centre points of ROIs for labels...');
for j = 1:length(i)
    vox    = find(roi==i(j));
    verts{j}  = v(vox,:);
    C(j,:) = spherefit(verts{j});
end
fprintf('  ... done! \n');
% % Plot the first roi, mark centre and label:
% scatter3(v(:,1),v(:,2),v(:,3),'k'); hold on
% scatter3(verts(:,1),verts(:,2),verts(:,3),'r')
% scatter3(C(:,1),C(:,2),C(:,3),'b*')

end