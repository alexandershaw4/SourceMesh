function [pos,overlay] = psudoAAL()

% Make psudo overlay from AAL

load AAL_labels.mat
load New_AALROI_6mm.mat

[SELECTION,OK] = listdlg('ListString',AAL_Labels);

Vals    = randi([-7 7],length(SELECTION),1);
overlay = zeros(5061,1);
pos     = template_sourcemodel.pos;

if OK
    for i = 1:length(SELECTION)
        this = find(all_roi_tissueindex==SELECTION(i));
        overlay(this) = Vals(i);
    end
end