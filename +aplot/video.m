function data = video(data,L,colbar,fpath,tv)
%

% OPTIONS
%--------------------------------------------------------------------------
num         = 1;   % number of brains, 1 or 2
interpl     = 1;   % interpolate
brainview   = 'T'; % [T]op, [L]eft or [R]ight
videolength = 10;  % length in seconds
extendvideo = 0;   % smooth/extend video by factor of

pos  = data.sourcemodel.pos;
mesh = data.mesh; 

data.video.opt.num         = num;
data.video.opt.interpl     = interpl;
data.video.opt.brainview   = brainview;
data.video.opt.videolength = videolength;
data.video.opt.extendvideo = extendvideo;

% Extend and temporally smooth video by linear interp between points
%--------------------------------------------------------------------------
if extendvideo > 0
    fprintf('Extending and smoothing video sequence by linear interpolation\n');
    time  = tv;
    for i = 1:size(L,1)
        dL(i,:) = interp(L(i,:),4);
    end
    L  = dL;
    tv = linspace(time(1),time(end),size(L,2));
end

data.video.t = tv;
%data.video.data = L;

% Overlay
%--------------------------------------------------------------------------
v  = pos;
x  = v(:,1);                    % AAL x verts
mv = mesh.vertices;             % brain mesh vertices
nv = length(mv);                % number of brain vertices
ntime = size(L,2);
try
    OL = zeros(size(L,1),nv,ntime); % this will be overlay matrix we average
catch
    fprintf('------------------------ ERROR ------------------------\n');
    fprintf('_______________________________________________________\n');
    fprintf('Projection matrix too big: M(Sources*MeshVertices*Time)\n');
    fprintf('M = size( %d , %d, %d )\n',size(L,1),nv,ntime);
    fprintf('Try: (1). Reducing source by using AAL template flag: ''template'',''aal'' \n');
    fprintf('     (2). Subsample time\n');
    return;
end
r  = 1200;                      % radius - number of closest points on mesh
r  = (nv/length(pos))*1.3;
w  = linspace(.1,1,r);          % weights for closest points
w  = fliplr(w);                 % 
M  = zeros( length(x), nv);     % weights matrix: size(len(mesh),len(AAL))
S  = [min(L)',max(L)'];

% Get centre point of cortical mesh so we know left/right
cnt = spherefit(mv);
Lft = mv(:,1) < cnt(1);
Rht = mv(:,1) > cnt(1);

Lft = find(Lft);
Rht = find(Rht);

% find closest points (assume both in mm)
%--------------------------------------------------------------------------
fprintf('Determining closest points between sourcemodel & template vertices\n');
for i = 1:length(x)

    % reporting
    if i > 1; fprintf(repmat('\b',[size(str)])); end
    str = sprintf('%d/%d',i,(length(x)));
    fprintf(str);    
    
    LR     = v(i,1);
    IsLeft = (LR-cnt(1)) < 0;
    
    if IsLeft; lri = Lft;
    else;      lri = Rht;
    end
    
    dist       = cdist(mv(lri,:),v(i,:));    
    [junk,ind] = maxpoints(dist,max(r,1),'min');   
    ind        = lri(ind);    
    OL(i,ind,:)= w'*L(i,:);
    M (i,ind)  = w;  
    
end
fprintf('\n');


if ~interpl
    OL = mean((OL),1); % mean value of a given vertex
else
    fprintf('Averaging local & overlapping vertices (wait...)');
    for i = 1:size(OL,2)
        for j = 1:size(OL,3)
            % average overlapping voxels 
            L(i,j) = sum( OL(:,i,j) ) / length(find(OL(:,i,j))) ;
        end
    end
    fprintf(' ...Done\n');
    OL = L;
end

% normalise and rescale
for i = 1:size(OL,2)
    this    = OL(:,i);
    y(:,i)  = S(i,1) + ((S(i,2)-S(i,1))).*(this - min(this))./(max(this) - min(this));
end

y(isnan(y)) = 0;
y  = full(y);

% spm mesh smoothing
fprintf('Smoothing overlay...\n');
for i = 1:ntime
    y(:,i) = spm_mesh_smooth(mesh, double(y(:,i)), 4);
end

data.video.data = y;
data.video.weights = M;

% close image so can reopen with subplots
if num == 2;
    close
    f  = figure;
    set(f, 'Position', [100, 100, 2000, 1000])
    h1 = subplot(121);
    h2 = subplot(122);
else
    switch brainview
        case 'T'; bigimg;view(0,90);
        case 'R'; bigimg;view(90,0);  
        case 'L'; bigimg;view(270,0); 
    end
    f = gcf;
end


% only project requested hemisphere
% switch data.hemi
%     case{'left','L','l'}; vi = data.mesh.vleft;
%     case{'right','R','r'};vi = data.mesh.vright;
%     otherwise;            vi = 1:length(data.mesh.vertices);
% end


% MAKE THE GRAPH / VIDEO
%--------------------------------------------------------------------------
try    vidObj   = VideoWriter(fpath,'MPEG-4');          % CHANGE PROFILE
catch  vidObj   = VideoWriter(fpath,'Motion JPEG AVI');
end

set(vidObj,'Quality',100);
set(vidObj,'FrameRate',size(y,2)/(videolength));
open(vidObj);

for i = 1:ntime
    
    if i > 1; fprintf(repmat('\b',[1 length(str)])); end
    str = sprintf('building: %d of %d\n',i,ntime);
    fprintf(str);
    
    switch num
        case 2
            plot(h1,gifti(mesh));
            hh       = get(h1,'children');
            set(hh(end),'FaceVertexCData',y(:,i), 'FaceColor','interp');    
            shading interp
            view(270,0);
            caxis([min(S(:,1)) max(S(:,2))]);
            material dull
            camlight left 

            plot(h2,gifti(mesh));
            hh       = get(h2,'children');
            set(hh(3),'FaceVertexCData',y(:,i), 'FaceColor','interp');    
            shading interp
            view(90,0);
            caxis([min(S(:,1)) max(S(:,2))]);
            material dull
            camlight right 
        
        case 1
            hh = get(gca,'children');
            set(hh(end),'FaceVertexCData',y(:,i), 'FaceColor','interp');
            caxis([min(S(:,1)) max(S(:,2))]);
            shading interp
    end
    
    try
        tt = title(num2str(tv(i)),'fontsize',20);
        P = get(tt,'Position') ;
        P = P/max(P(:));
        set(tt,'Position',[P(1) P(2)+70 P(3)]) ;
    end
    
    set(findall(gca, 'type', 'text'), 'visible', 'on');
    
    if colbar
        drawnow; pause(.5);
        a1 = gca;
        axb = axes('position', get(a1, 'position'));
        set(axb,'visible','off')
        axes(axb);
        colorbar('peer',a1,'South');
    end
    drawnow;
            
              

    currFrame = getframe(f);
    writeVideo(vidObj,currFrame);
end
close(vidObj);


    
end