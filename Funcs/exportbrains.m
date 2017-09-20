function exportbrains(fname)
% exportbrains

view(0,90);  drawnow; print([fname '_top']   ,'-dpng','-r600') ;
view(270,0); drawnow; print([fname '_left']  ,'-dpng','-r600') ;
view(90,0);  drawnow; print([fname '_right'] ,'-dpng','-r600') ;
view(0,0);   drawnow; print([fname '_back']  ,'-dpng','-r600') ;
view(0,-90); drawnow; print([fname '_bottom'],'-dpng','-r600') ;
view(0,-90); drawnow; print([fname '_bottom'],'-dpng','-r600') ;
view(-180,0);drawnow; print([fname '_front'] ,'-dpng','-r600') ;