function data = addlabels(data,V,all_roi_tissueindex,thelabels)
% Add labels to the plot.
%
% If using AAL90 sourcemodle, these are automatic.
%
% If using another sourcemodel:
% - provide the all_roi_tissueindex from fieldtirp. This is a
% 1xnum_vertices vector containing indices of rois (i,e. which verts belong
% to which rois).
% Also provide labels!
%
pos = data.sourcemodel.pos;

if ( ~isempty(thelabels) && ~isempty(all_roi_tissueindex) ) &&...
   ( length(pos) == length(all_roi_tissueindex) ) &&...
   ( length(thelabels) == length(unique(all_roi_tissueindex(all_roi_tissueindex~=0))) )
    
    labels = strrep(thelabels,'_',' ');
    v      = aplot.get_roi_centres(pos,all_roi_tissueindex);
    roi    = all_roi_tissueindex;
    
elseif length(V) == 90
    
    load('AAL_labels');
    labels = strrep(labels,'_',' ');
    v      = pos*0.95;
    roi    = 1:90;
elseif (length(V) == length(thelabels)) &&...
       (length(V) == length(pos))
    
   labels = strrep(thelabels,'_',' ');
    v = pos*0.95;
    roi = 1:length(V);
else
    fprintf('Labels info not right!\n');
    return
end

data.labels.roi     = roi;
data.labels.labels  = labels;
data.labels.centres = v;

% compile list of in-use node indices
%--------------------------------------------------------------------------
to = []; from = []; V(isnan(V)) = 0;
for i  = 1:size(V,1)
    ni = find(logical(V(i,:)));
    if any(ni)
        to   = [to   roi(ni)];
        from = [from roi(repmat(i,[1,length(ni)])) ];
    end
end

AN  = unique([to,from]);
AN  = AN(AN~=0);
off = 1.5;
data.labels.in_use = AN;

if sum(v(:,3)) == 0
     off3 = 0;
else;off3 = off;
end

% add these to plot with offset
%--------------------------------------------------------------------------
for i = 1:length(AN)
    L = labels{AN(i)};
    switch L(end)
        case 'L';
            t(i) = text(v(AN(i),1)-(off*5),v(AN(i),2)-(off*5),v(AN(i),3)+off3,L);
        case 'R';
            t(i) = text(v(AN(i),1)+(off*2),+v(AN(i),2)+(off*2),v(AN(i),3)+off3,L);
        otherwise
            t(i) = text(v(AN(i),1),v(AN(i),2),v(AN(i),3),L);
    end
end
try
    % this fails if the network was empty!
    set(t,'Fontsize',14)
end

end