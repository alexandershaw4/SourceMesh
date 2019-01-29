function [data,i] = sort_template(data,i)
% If specified a template model, register data to it and return splined data as
% well as weights

if ~isfield(i,'pos')
    i.pos = data.sourcemodel.pos;
end
try
    data.template.model  = i.model;
    data.template.labels = i.labels;
end
if i.template
    atlas = dotemplate(i.model);
    rois  = get_roi_centres(atlas.template_sourcemodel.pos,atlas.all_roi_tissueindex);
    
    atlas.template_sourcemodel.pos = rois;
    atlas = rmfield(atlas,'all_roi_tissueindex');
    
    reg = interp_template(data.sourcemodel,rois);
    atlas.M    = reg.M;
    data.atlas = atlas;
    NM         = atlas.M;
    
    % rescale so not change amplitudes
    m  = max(NM(:));
    NM = NM/m; 
    
    % update sourcemodel and labels
    data.sourcemodel = atlas.template_sourcemodel;
    if i.labels; i.thelabels = atlas.AAL_Labels; end
    
    % overlay data
    if isfield(i,'L')
        if isnumeric(i.L) && ndims(i.L) ~= 3
            S  = [min(i.L(:)) max(i.L(:))];
            NL = i.L(:)'*NM;
            L  = S(1) + ((S(2)-S(1))).*(NL - min(NL))./(max(NL) - min(NL));
            L(isnan(L))=0;
            i.L = L;
        end
    end
    
    % network
    if isfield(i,'A')
        if isnumeric(i.A)
            S  = [min(i.A(:)) max(i.A(:))];
            NL = NM'*i.A*NM;
            A  = S(1) + ((S(2)-S(1))).*(NL - min(NL(:)))./(max(NL(:)) - min(NL(:)));
            A(isnan(A)) = 0;
            i.A = A;
        end
    end
    
    % video data
    if isfield(i,'V')
        S  = [min(i.V(:)) max(i.V(:))];
        for j = 1:size(i.V,2) % over time points
            NL(:,j) = i.V(:,j)'*NM;
        end
        V  = S(1) + ((S(2)-S(1))).*(NL - min(NL))./(max(NL) - min(NL));
        V(isnan(V))=0;
        if orthog
            % dont use this
            V = symm_orthog(V);
        end
        V(isnan(V))=0;
        i.V = V;
    end
        
end

end