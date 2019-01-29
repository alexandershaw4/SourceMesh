function atlas = dotemplate(model)
% Put dense sourcemodel into an atlas space using ICP and linear
% interpolation
%
%
%

switch model
    case lower({'aal','aal90'});   load New_AALROI_6mm.mat
    case lower('aal58');           load New_58cortical_AALROI_6mm
    case lower('aal78');           load New_AALROI_Cortical78_6mm
    otherwise
        fprintf('Model not found.\n');
        return;
end

atlas.AAL_Labels = AAL_Labels;
atlas.all_roi_tissueindex = all_roi_tissueindex;
atlas.template_sourcemodel = template_sourcemodel;

end