function data = parse_labels(i,data)
% decide which labels to include depending on what we're plotting
if i.labels 
    if     isfield(i,'A')
                if isnumeric(i.A)
                    data = addlabels(data,i.A,i.all_roi_tissueindex,i.thelabels);
                elseif ischar(i.A)
                    E = data.network.edge;
                    data = addlabels(data,E,i.all_roi_tissueindex,i.thelabels);
                end
           
    elseif isfield(i,'N')
        if sum(ismember(size(i.N),[1 90])) == 2
            data = addlabels(data, diag(i.N),i.all_roi_tissueindex,i.thelabels);
        elseif sum(ismember(size(i.N),[1 90])) == 1
            data = addlabels(data, diag(sum(i.N,2)),i.all_roi_tissueindex,i.thelabels);
        end
        
    else;  n    = length(data.sourcemodel.pos);
           data = addlabels(data, ones(n,n),i.all_roi_tissueindex,i.thelabels);
    end
end
end