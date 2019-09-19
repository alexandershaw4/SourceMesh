function im2vid(fn,varargin)
% Capture (current) 3D figure, make video by rotating camera angle 360
% degrees. 
%
% Usage: im2vid('myvideo')
%
% Other useage:
%        im2vid('myvideo','quick');
%        im2vid('myvideo','background');
%
% Quick mode is fast, but requires that the figure window be showing and 
% on top of all other windows. 
% Background is slower but you can carry on doing other stuff - it writes
% each frame to a temporary png.
%
% AS

if nargin > 1;
    method = varargin{1};
else
    method = 'background'; % 'quick' or 'background'
end

[fp,fn,fe] = fileparts(fn);
if isempty(fp)
    try   fpath = [pwd, '/', fn];
    catch fpath = [pwd, '/im'];
    end
else
    fpath = [fp '/' fn];
end

pathname = pwd;
D        = 'right';
cam      = camlight(D);
vidObj   = VideoWriter(fpath,'MPEG-4');
set(vidObj,'Quality',100);
open(vidObj);

for num = 1:360
    if num > 1; fprintf(repmat('\b',[1 length(str)])); end
    str = sprintf('building: %d of %d\n',num,360);
    fprintf(str);
    camorbit(1,0,'camera');
    camlight(cam); drawnow;
    
    switch method
        case 'quick';
            currFrame = getframe(gcf);
            writeVideo(vidObj,currFrame);
        case 'background'
            
            print(gcf,[pathname,'/temp.png'],'-dpng','-r600');
            %export_fig([pathname,'/temp.png'],'-transparent','-m2','-nocrop');
            tempimg   = imread([pathname,'/temp.png']);

            currFrame = im2frame(tempimg);
            currFrame.cdata = imresize(currFrame.cdata,[1068,1470]);
            writeVideo(vidObj,currFrame);
    end
end

close(vidObj);
switch method
    case 'background'
        delete([pathname,'temp.png']);
end
fprintf('finished\n');