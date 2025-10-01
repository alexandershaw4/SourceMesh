
% The spm file that has been souce localised (e.g. has D.inv{1}):
D = spm_eeg_load('bawinica_ef1_30spmeeg_faces_PSP_0008_1_sss')

% The forward model mesh (D.inv{1}.forward)
mesh = D.inv{1}.forward(1).mesh;
v    = mesh.vert;

% Extract inverse components for reconstruction:
Mods = D.inv{1}.inverse.modality;
U    = D.inv{1}.inverse.U;
Ic   = D.inv{1}.inverse.Ic;
S    = D.inv{1}.inverse.scale;
M    = D.inv{1}.inverse.M;

% Average trials, D0 = (Chans * time)
D0   = mean(D(:,:,:),3);

for i  = 1:length(Mods)
    fprintf('Reconstructing modality: %s\n',Mods{i});
  % R    =  (spatial modes * wights) * Channel Data
    R{i} = ( U{i}          * S(i)  ) * D0(Ic{i},:);
end

% S = Filters * R(spatially reduced channel data by time)
Source = M*cat(1,R{:});
time   = D.time;

% Plot the average activity on the forward model mesh:
atemplate('gifti',gifti(mesh),'overlay',mean(Source,2))

% Align & plot the average activity on a smooth mesh:
atemplate('sourcemodel',v,'overlay',mean(Source,2))

% Plot video of activity on the forward model mesh, ...
% but put in AAL90 space for faster computation:
atemplate('gifti',gifti(mesh),'sourcemodel',v,'template','aal90','video',Source,'myvid',time)

% Plot video of activity on the smooth mesh, ...
% but put in AAL90 space for faster computation & subsample time:
atemplate('sourcemodel',v,'template','aal90','video',Source(:,1:6:end),'myvid2',time(1:6:end))
