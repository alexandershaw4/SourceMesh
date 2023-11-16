function [pos,labs] = map_pos_to_aal(pos,vi)


aal = load('SuperAAL.mat');
labels = load('AAL_labels.mat');

% alginmane and scaling
pos = alignclouds(pos,aal.v);

% now find centre of each parcel in input data
roi = unique(vi);

for i = 1:length(roi)
    I = find(vi == roi(i));
    
    verts = pos(I,:);

    C = spherefit(verts);

    if any(isnan(C)) || any(isinf(C))
        C = mean(verts,1);
    end

    npos(i,:) = C;
end

%npos = npos - (pinv(cdist(npos,npos))*npos);

for i = 1:size(npos,1)

    % closest in euclidean space
    [~,ind] = min(cdist(npos(i,:),aal.v));

    labs{i} = labels.labels{aal.vi(ind)};


end

% output
L = {...
    'Cuneus_L'
    'Cuneus_R'
    'Occipital_Inf_L'
    'Occipital_Inf_R'
    'Postcentral_L'
    'Postcentral_R'
    'Rolandic_Oper_L'
    'Rolandic_Oper_R'
    'Precentral_L'
    'Precentral_R'
    'Parietal_Sup_L'
    'Parietal_Sup_R'
    'Temporal_Sup_L' %13
    'Temporal_Sup_R' %14
    'Angular_L'
    'Angular_R'
    'Temporal_Pole_Sup_L'
    'Temporal_Pole_Sup_R'
    'Paracentral_Lobule_L'
    'Paracentral_Lobule_R'
    'Parietal_Inf_L'
    'Parietal_Inf_R'
    'Frontal_Inf_Tri_L'
    'Frontal_Inf_Tri_R'
    'Occipital_Mid_L'
    'Occipital_Mid_R'
    'Frontal_Mid_L'
    'Frontal_Mid_R'
    'Frontal_Sup_L' % 29 yes sup
    'Frontal_Sup_R' % 30
    'Frontal_Mid_Orb_L'
    'Frontal_Mid_Orb_R'
    'Temporal_Mid_L' % 33
    'Temporal_Mid_R' % 34
    'Frontal_Sup_Orb_L' % 35 
    'Frontal_Sup_Orb_R' % 36
    'Cingulum_Post_L'
    'Cingulum_Ant_L'};