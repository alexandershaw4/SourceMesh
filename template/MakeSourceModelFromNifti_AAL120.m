

%mri = ft_read_mri('/Users/Alex/Dropbox/AAL_Template/AAL_Labels.nii');
mri = ft_read_atlas('~/fieldtrip-master/template/atlas/aal/ROI_MNI_V4.nii');
ana = mri.tissue;

S = size(ana);
x = 1:S(1); 
y = 1:S(2); 
z = 1:S(3); 

% Map volume to a Cartesian grid with indicies for vertex membership
v  = [];
vi = [];

ID = unique(ana(:));

for i = 1:length(ID)  %90
    if ID(i) > 0
        If = find(ana==ID(i));
        [xi,yi,zi] = ind2sub(S,If);

        newv = [x(xi); y(yi); z(zi)]';

        if any(If(:))

            % assess proximity of verts in this group
            t = delaunay(newv); t = t(:,1:3);
            [nf,nv] = reducepatch(maker(t,newv),0.1);
            newv = nv;
            newi = ones(length(newv),1)*i;

            v  = [v;   newv];
            vi = [vi ; newi];
        end
    end
end

c = spherefit(v);
v = v - repmat(c,[size(v,1),1]);

%save('~/code/MeshAAL/template/DenseAAL','v','vi')
%save('~/code/MeshAAL/template/HarvOx','v','vi')
labels = mri.tissuelabel';
save('~/code/MeshAAL/template/AAL120','v','vi','labels')

