
PlotStat = 'Vis100B_Av.nii'; % The nifti overlay file

afigure;   % Opens a large figure window

D = atemplate('mesh','def4','inflate','overlay',{'curvature' PlotStat},...
    'method','raycast','thresh',6,'hemi','l');

axis equal;

aJointColorbar(D); % Adds a colorbar BUT, CAXIS WONT WORK... :(