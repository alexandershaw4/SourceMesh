% Example Rendering Overlays.....

load New_AALROI_6mm.mat   % sourcemodel
load SaveVolumeTstats.mat % file containing 'Plotstat': vector length(5061)

projmeth = 'raycast';

PlotStat = squeeze(sum(sum(Tstat,1),2));

% Just render PlotStat as colours onto the brain (no curvature)
afigure;
D = atemplate('mesh','def5','overlay', PlotStat,...
            'sourcemodel',template_sourcemodel.pos,...
            'method',projmeth);
 
axis equal;
caxis([-3 3]);
colormap(cmocean('balance'));

% Check it against a 3d scatter with overlay of the original positions:
pp = template_sourcemodel.pos;
figure,scatter3(pp(:,1),pp(:,2),pp(:,3),90,PlotStat,'filled');
view([0 90]);


% OR WITH CURVATURE & COLORBAR ALSO, USING EXPLICIT THRESH:
%------------------------------------------------------------
afigure;
D = atemplate('mesh','def4','inflate','overlay',{'curvature' PlotStat},'sourcemodel',template_sourcemodel.pos,...
    'method',projmeth,'thresh',3);axis equal;

aJointColorbar(D); % BUT, CAXIS WONT WORK... :(
