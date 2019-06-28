% Make a parcellation that uses Harvard-Oxford Cortical & AAL Cerebellar regions
%-------------------------------------------------------------------------------

AAL = load('AAL116.mat');
HOA = load('HarvOx.mat');

% coregister both points to MNI models so they're approx aligned:
% register HarvOx to Cortex2
% register AAL (with Cereb) to Cortex3 (with Cereb), both MNI
figure;D = atemplate('mesh','def2','sourcemodel',HOA.v);
HOA.v = D.sourcemodel.pos; close;

figure;D = atemplate('mesh','def3','sourcemodel',AAL.v);
AAL.v = D.sourcemodel.pos; close;



Cereb = 91:116;
CLab  = AAL.labels(Cereb);
C     = [];
VI    = [];

for j  = 1:length(Cereb)
    i  = find( AAL.vi == Cereb(j) );
    C  = [C  ; AAL.v(i,:) ];
    VI = [VI ; (i*0)+j ];
end


% plot conbination of models on mesh:
afigure;D = atemplate('mesh','def3');

% Harvard-Oxford cortex
scatter3(HOA.v(:,1),HOA.v(:,2),HOA.v(:,3),'b','filled'); hold on;

% AAL Cerebellum & Vermis
scatter3(C(:,1),C(:,2),C(:,3),'r','filled');

% Combine model vertices and optimise alignment with MNI
mixmodel    = [HOA.v  ; C ];
mixmodel_iv = [HOA.vi ; VI + max(HOA.vi)];
mixmodel_al = align_clouds_3d(D.mesh.vertices,mixmodel);

% Now move each point to the closest mesh point
D0 = cdist(mixmodel_al,D.mesh.vertices);

for i = 1:size(D0,1)
    [~,p] = min(D0(i,:));
    mixmodel_al_re(i,:) = D.mesh.vertices(p,:);
end
    

close all;

v  = mixmodel_al_re;
vi = mixmodel_iv; 

% also keep a record of which vertices came from which model:
v_source( 1:length(HOA.vi))                 = 0;
v_source( length(HOA.vi) + (1:length(VI)) ) = 1;

Cort = v(v_source==0,:);
Cere = v(v_source==1,:);

% re-render aligned models:
afigure;D = atemplate('mesh','def3'); hold on;

scatter3(Cort(:,1),Cort(:,2),Cort(:,3),50,'b','filled');
scatter3(Cere(:,1),Cere(:,2),Cere(:,3),50,'r','filled');

% save it
labels = [HOA.labels(:,2) ; CLab ];
save('/Users/Alex/code/MeshAAL/template/HOA_Cerebellum','v','vi','Cort','Cere','labels','v_source');

