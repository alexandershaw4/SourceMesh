

%mri = ft_read_mri('/Users/Alex/Dropbox/AAL_Template/AAL_Labels.nii');
%mri = ft_read_mri('HarvardOxford-cort-maxprob-thr50-1mm.nii');
mri = ft_read_mri('/Users/Alex/fieldtrip-master/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii');
ana = mri.anatomy;

S = size(ana);
x = 1:S(1); 
y = 1:S(2); 
z = 1:S(3); 

% Map volume to a Cartesian grid with indicies for vertex membership
v  = [];
vi = [];
VALS = unique(ana(:));
VALS(VALS==0)=[];
for i = 1:length(VALS) %48  %90
    If = find(ana==i);
    [xi,yi,zi] = ind2sub(S,If);
    
    newv = [x(xi); y(yi); z(zi)]';
    
    if any(If(:))
    
        % assess proximity of verts in this group
        t = delaunay(newv); t = t(:,1:3);
        [nf,nv] = reducepatch(maker(t,newv),0.2);
        %nf=t;nv=newv;
        newv = nv;
        newi = ones(length(newv),1)*i;

        v  = [v;   newv];
        vi = [vi ; newi];
    end
end

c = spherefit(v);
v = v - repmat(c,[size(v,1),1]);

% optionally compute face and recude using points_to_alphaShape
[f,v,Ci] = points_to_alphashape(v);
vi = vi(Ci);
c = spherefit(v);
v = v - repmat(c,[size(v,1),1]);

% separate hemispheres if labelled the same in each hemisphere - e.g. the
% HCP360/250 atlasses

ui = unique(vi);
for i = 1:length(ui)

    x = find(vi==ui(i));

    L = find(v(x,1)<0);
    R = find(v(x,1)>0);
    
    % rename
    vi(x(L)) = i;
    vi(x(R)) = i + 180;
end


%save('~/code/MeshAAL/template/DenseAAL','v','vi')
save('~/Dropbox/code/MeshAAL/template/HCP360','v','vi','f')

% trisurf(f,v(:,1),v(:,2),v(:,3),'facevertexcdata',vi);
% colormap(alexcmap); axis equal

% lighter or sparser-version
%-------------------------------------------
mri = ft_read_mri('/Users/Alex/Dropbox/AAL_Template/AAL_Labels.nii');
ana = mri.anatomy;

S = size(ana);
x = 1:S(1); 
y = 1:S(2); 
z = 1:S(3); 

% Map volume to a Cartesian grid with indicies for vertex membership
v  = [];
vi = [];
for i = 1:90
    If = find(ana==i);
    [xi,yi,zi] = ind2sub(S,If);
    
    newv = [x(xi); y(yi); z(zi)]';
    
    % assess proximity of verts in this group
    t = delaunay(newv); t = t(:,1:3);
    [nf,nv] = reducepatch(maker(t,newv),0.05);
    newv = nv;
    
    newi = ones(length(newv),1)*i;
    
    v  = [v;   newv];
    vi = [vi ; newi];
end

c = spherefit(v);
v = v - repmat(c,[size(v,1),1]);

save('~/code/MeshAAL/template/LightAAL','v','vi')