function [pos,overlay] = psudoAAL(type)
% Make psudo overlay from AAL90 sourcemodel
%
% 1) call [pos,overlay] = psudoAAL('r') to select nodes to be non-zero from a
% popup list. 'r' indicates the values for these nodes will be randomly
% generated between -7 - 7.
%
% 2) call [pos,overlay] = psudoAAL()  with no input to select nodes to be 
% non-zero from a popup list. No input indicates values assigned to
% selected nodes will increase from 1:num selected.
%
% AS


load AAL_labels.mat
load New_AALROI_6mm.mat

if nargin < 1
    type = ' ';
end


[SELECTION,OK] = listdlg('ListString',AAL_Labels);

switch type
    case 'r';
        Vals    = randi([-7 7],length(SELECTION),1);
    otherwise
        Vals    = 1:length(unique(SELECTION));
end
    
overlay = zeros(5061,1);
pos     = template_sourcemodel.pos;

if OK
    for i = 1:length(SELECTION)
        this = find(all_roi_tissueindex==SELECTION(i));
        overlay(this) = Vals(i);
    end
end