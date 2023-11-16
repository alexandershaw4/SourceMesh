
PlotStat = 'Vis100B_Av.nii'; % The nifti overlay file

afigure;   % Opens a large figure window

D = atemplate('mesh','def5','inflate','overlay',{'curvature' PlotStat},...
    'method','raycast','thresh',4,'hemi','l');

axis equal;

aJointColorbar(D); % Adds a colorbar BUT, CAXIS WONT WORK... :(