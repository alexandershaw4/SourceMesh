
% PlotStat = 'Vis100B_Av.nii'; % The nifti overlay file
% 
% afigure;   % Opens a large figure window
% 
% % D = atemplate_patched('mesh','def5','inflate','overlay',{'curvature' PlotStat},...
% %     'method','raycast','thresh',1);
% 
% data.overlay_opts = struct('normal_depth_mm',2.5,'normal_steps',7,'aggregate','max', ...
%                            'inpaint_iters',1,'smooth_iters',3);
% 
% 
% D = atemplate_patched('mesh','def5','inflate','overlay',PlotStat,...
%     'method','raycast');
% 
% 
% 
% axis equal;
% 
%aJointColorbar(D); % Adds a colorbar BUT, CAXIS WONT WORK... :(


% New - much better raycast alignment

PlotStat = 'Vis100B_Av.nii';
%PlotStat = 'nTEST_Gamma,60-90,60-90Hz_pairedt.nii';
PlotStat = 'nTEST_Beta_desynch,15-30,15-30Hz_pairedt.nii';

opts = struct('normal_depth_mm',2.5,'normal_steps',7,'aggregate','max', ...
              'inpaint_iters',1,'smooth_iters',3);

afigure;
D = atemplate('mesh','def5','inflate', ...
              'overlay',PlotStat, ...
              'overlay_opts',opts, ...
              'method','raycast');
axis equal;

