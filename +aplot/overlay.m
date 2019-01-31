function data = overlay(data,L,write,fname,colbar)
% Functional overlay plotter
%
% mesh is the gifti / patch
% L is the overlay (90,1)
% write is boolean flag
% fname is filename is write = 1;
%


% Add this special case, where using default 81k mesh and 90-node AAL
% overlay, we'll use pre-computed weights for speed
%--------------------------------------------------------------------------
if isnumeric(L) && ndims(L)==2 && length(L)==90 && length(data.mesh.vertices)== 81924 && ...
        ischar(data.overlay.method) && strcmp(data.overlay.method,'euclidean')
    
     fprintf('Using default AAL90 weights for this mesh\n');
     load('AAL90DefaultWeights','M','NumComp','indz','w');
     
     % incorporate overlay into precomputed weights matrix
     for k = 1:length(L)
         M(k,indz(k,:)) = L(k) * M(k,indz(k,:));
     end
     
     % normalise by number of overlapping points at this vertex
     for i = 1:size(M,2)
         y(i) = sum( M(:,i) ) / length(find(M(:,i))) ;
     end   
     
     % rescale y by L limits
     S  = [min(L(:)),max(L(:))];   
     y  = S(1) + ((S(2)-S(1))).*(y - min(y))./(max(y) - min(y));
     y(isnan(y)) = 0;
     L  = y;
     data.method = 'precomputed (AAL)';
end


% deal with the overlay if it's a filename or volume
%-------------------------------------------------------------
if ~isnumeric(L) || (isnumeric(L) && ndims(L)==3)
    % is this is filename of a nifti or gifti file
   [L,data] = aplot.parse_overlay(L,data);
   
   if isempty(L)
        fprintf('Overlay does not match sourcemodel!\n');
        return;
   end
   if isfield(data,'template')
       if isfield(data.template,'model')
           fprintf('Doing atlas registration\n');
           i.template = 1;
           i.model    = data.template.model;
           i.labels   = data.template.labels;
           i.L        = L;
           [data,i]   = aplot.sort_template(data,i);
           L          = i.L;
       end
   end
end

% method for searching between the 3D coordinate systems
%-------------------------------------------------------------
if ischar(data.overlay.method)
    if ismember(lower(data.overlay.method),{'euclidean','spheres','precomputed (AAL)','raycast','aal','aal_light'})
         method = data.overlay.method;
    else,method = 'euclidean';  
    end
else
    method = data.overlay.method{1};
end

% If atlas data and peaks frequested, label them
%-------------------------------------------------------------
data.overlay.orig = L;
if data.overlay.peaks
    n     = mean(L)+(2*std(L));
    [V,I] = find(abs(L) > n);

    if isfield(data,'atlas')
        Lab = data.atlas.AAL_Labels;
        data.overlay.Peaks.Labels = Lab(I);
        data.overlay.Peaks.Values = L(I);
    end
    
end

% interp shading between nodes or just use mean value?
%--------------------------------------------------------------------------
interpl = 1; 
pos     = data.sourcemodel.pos;
mesh    = data.mesh;



% if overlay,L, is same length as mesh verts, just plot!
%--------------------------------------------------------------------------
if length(L) == length(mesh.vertices)
    fprintf('Overlay already fits mesh! Plotting...\n');
    
    % spm mesh smoothing
    fprintf('Smoothing overlay...\n');
    y = spm_mesh_smooth(mesh, double(L(:)), 4);
    percNaN = length(find(isnan(L)))/length(L)*100;
    newpNaN = length(find(isnan(y)))/length(y)*100;
    
    % when using a NaN-masked overlay, smoothing can result in all(nan) or
    % an increase in the number of nans: enforce a 5% tolerance on this, which
    % forces reverting to the uns-smoothed version if reached
    if all(isnan(y)) || newpNaN > (percNaN*1.05)
        fprintf('Reverting to non-smoothed overlay due to too many NaNs\n');
        y = L(:);
    end
    
    set(mesh.h,'FaceVertexCData',y(:),'FaceColor','interp');
    
    drawnow;
    shading interp
    % force symmetric caxis bounds
    s = max(abs(y(:))); caxis([-s s]);
    colormap('jet');
    alpha 1;
    
    if colbar
        %colorbar('peer',a1,'South');
        data.overlay.cb = InteractiveColorbar;
    end
    
    if write == 1
        fprintf('Writing overlay gifti file: %s\n',[fname 'Overlay.gii']);
        g       = gifti;
        g.cdata = double(y);
        g.private.metadata(1).name  = 'SurfaceID';
        g.private.metadata(1).value = [fname 'Overlay.gii'];
        save(g, [fname  'Overlay.gii']);
    elseif write == 2
            fprintf('Writing mesh and overlay as STL object\n');
        % write STL
        m.vertices = double(mesh.vertices);
        m.faces    = double(mesh.faces);
        y          = double(y);
        cdata      = mean(y(m.faces),2);

        % Write binary STL with coloured faces
        cLims = [min(cdata) max(cdata)];      % Transform height values
        nCols = 255;  cMap = jet(nCols);    % onto an 8-bit colour map
        fColsDbl = interp1(linspace(cLims(1),cLims(2),nCols),cMap,cdata);

        fCols8bit = fColsDbl*255; % Pass cols in 8bit (0-255) RGB triplets
        stlwrite([fname '.stl'],m,'FaceColor',fCols8bit)
        
    elseif write == 3 
        % write vrml
        fprintf('Writing vrml (.wrl) 3D object\n');
        vrml(gcf,[fname]);
        
    elseif write == 4
%         fprintf('Generating nifti volume for writing\n');
%         V = sm2vol(mesh.vertices,256,y,50);
%         fprintf('Writing nifti volume\n');
%         niftiwrite(V, [fname '.nii']);
    end


else

% otherwise find closest points (assume both in mm)
%--------------------------------------------------------------------------

% Overlay
v  = pos;                       % sourcemodel vertices
x  = v(:,1);                    % AAL x verts
mv = mesh.vertices;             % brain mesh vertices
nv = length(mv);                % number of brain vertices
S  = [min(L(:)),max(L(:))];     % min max values

switch method
    case{'raycast'}
    otherwise
    if write == 2
        r = (nv/length(pos))*5;
        w  = linspace(.1,1,r);          % 
        w  = fliplr(w);                 % 
    elseif write == 3
        r = (nv/length(pos))*3;
        w  = linspace(.1,1,r);          % 
        w  = fliplr(w);                 % 
    end
end

% Get centre point of cortical mesh so we know left/right for hemispheric
% separation
cnt = spherefit(mv);
Lft = mv(:,1) < cnt(1);
Rht = mv(:,1) > cnt(1);
Lft = find(Lft);
Rht = find(Rht);

fprintf('Determining closest points between sourcemodel & template vertices\n');
mr = mean(mean(abs(mv)-repmat(spherefit(mv),[size(mv,1),1])));


% Switch which projection method to use:
%--------------------------------------------------------------------------
switch lower(method)
    
    case 'raycast'
        % This is an attempt at employing the ray casting method
        % implemented in mri3dX. It grids the functional data and searches
        % along each face's normal from -1.5 to 1.5 in 0.05 mm steps.
        %
        % The functional overlay vector returned in the overlay substructure 
        % has one value per FACE of the mesh, although this is converted to
        % vertex values by taking the maximum value of each triangle
        
        wb = waitbar(0,'Ray casting: Please wait...');
        
        % Ray cast from FACES or from VERTICES: SET 'face' / 'vertex'
        UseFaceVertex = 'face'; 
        RND = 1;
        
        % Grid resolution
        nmesh.vertices = data.mesh.vertices * .5;
        dv             = v * .5;
                
        % make new mesh and overlay points, decimated / rounded to integers (mm)
        nmesh.vertices = double(round(nmesh.vertices*RND)/RND);
        nmesh.faces    = double(data.mesh.faces);
        dv             = round(dv*RND)/RND;
           
        waitbar(.2,wb,'Ray casting: Gridding data');
        
        % volume the data so vertices are (offset) indices
        fprintf('Gridding data for ray cast\n');
        vol = zeros( (max(dv) - min(dv))+1 );
        ndv = min(dv)-1;
        
        for i = 1:length(dv)
            if L(i) ~= 0
                a(1)  = L(i);
                a(2)  = vol(dv(i,1)-ndv(1),dv(i,2)-ndv(2),dv(i,3)-ndv(3));
                [~,I] = max(abs(a));
                vol(dv(i,1)-ndv(1),dv(i,2)-ndv(2),dv(i,3)-ndv(3)) = a(I);                
            end
        end
                
        waitbar(.4,wb,'Ray casting: Smoothing');
        
        % Smooth volume
        fprintf('Volume Smoothing & Rescaling  ');tic        
        vol  = smooth3(vol,'box',3);        
        V    = spm_vec(vol);
        V    = S(1) + (S(2)-S(1)).*(V(:,1) - min(V(:,1)))./(max(V(:,1)) - min(V(:,1)));
        vol  = spm_unvec(V, vol); 
        fprintf('-- done (%d seconds)\n',round(toc)); 
        
        switch UseFaceVertex
            
            case 'face'
                
                % Load or compute FACE normals and centroids
                %----------------------------------------------------------
                if length(mv) == 81924
                    % use precomputed for deault mesh
                    load('DefaultMeshCentroidsNormals','FaceCent','FaceNorm')
                    fprintf('Using precomputed centroids & normals for default mesh\n');
                    f = nmesh.faces;
                else
                    
                    waitbar(.6,wb,'Ray casting: Computing face norms and centroids');
                    
                    % Compute face normals
                    %------------------------------------------------------
                    fprintf('Computing FACE Normals & Centroids  '); tic;
                    tr = triangulation(nmesh.faces,nmesh.vertices(:,1),...
                                        nmesh.vertices(:,2),nmesh.vertices(:,3));
                    FaceNorm = tr.faceNormal;

                    % Compute triangle centroids
                    %------------------------------------------------------
                    f        = nmesh.faces;
                    for If   = 1:length(f)
                        pnts = [nmesh.vertices(f(If,1),:);... 
                                nmesh.vertices(f(If,2),:);...
                                nmesh.vertices(f(If,3),:)];
                            
                        % Triangle centroid
                        FaceCent(If,:) = mean(pnts,1);
                    end
                    
                    fprintf('-- done (%d seconds)\n',round(toc));
                end
                
                % If a depth vector was specified use that, otherwise
                % deault
                if isfield(data.overlay,'depth') && ~isempty(data.overlay.depth)
                      step = data.overlay.depth;
                else; step   = -1.5:0.05:1.5;
                end
                fprintf('Using depths: %d to %d mm in increments %d\n',...
                    step(1), step(end), round((step(2)-step(1))*1000)/1000 );
                fcol   = zeros(length(step),length(f));
                
                
            case 'vertex'
                
                % Compute VERTEX normals
                fprintf('Computing VERTEX normals\n');
                FaceNorm = spm_mesh_normals(nmesh,1);
                
                % In this case, centroids are the vertices themselves
                FaceCent = nmesh.vertices;
                
                step    = -1.5:0.05:1.5;
                fcol    = zeros(length(step),length(mv));
        end
    
        % Now search outwards along normal line
        %-----------------------------------------------------------------
        waitbar(.8,wb,'Ray casting: casting');
        
        nhits  = 0; tic    ;
        perc   = round(linspace(1,length(step),10));
        for i  = 1:length(step)
            
            % keep count of num hits
            hits{i} = 0;
            
            % print progress
            if ismember(i,perc)
                fprintf('Ray casting: %d%% done\n',(10*find(i==perc)));
            end
            
            % the new points
            these = FaceCent + (step(i)*FaceNorm);
            
            % convert these points to indices of the volume
            these(:,1) = these(:,1) - ndv(1);
            these(:,2) = these(:,2) - ndv(2);
            these(:,3) = these(:,3) - ndv(3);
            these      = round(these*RND)/RND;
            
            % retain these
            FaceNormLine(i,:,:) = these;
            
            % values at volume indices
            for j = 1:length(these)
                try
                    fcol(i,j) = vol(these(j,1),these(j,2),these(j,3));
                    hits{i}   = hits{i} + 1;
                end
            end
        end

        fprintf('Finished in %d sec\n',round(toc));
        
        % Retain largest absolute value for each face (from each depth)
        [~,I] = max(abs(fcol));
        for i = 1:length(I)
            nfcol(i) = fcol(I(i),i);
        end
        fcol = nfcol;
        
        % add the values - either 1 per face or 1 per vertex - to the mesh
        %------------------------------------------------------------------
        waitbar(.9,wb,'Ray casting: sorting...');
        switch UseFaceVertex
            case 'face'
                % Set face colour data on mesh, requires setting FaceColor= 'flat'
                
                %set(mesh.h,'FaceVertexCData',fcol(:));
                %mesh.h.FaceColor = 'flat';
                
                % Or calculate vertex maxima from faces and use interp'd
                % vertex colouring
                f  = nmesh.faces;
                ev = mv*0;
                
                % these are the vals at the three corners of each triangle
                ev(f(:,1),1) = fcol;
                ev(f(:,2),2) = fcol;
                ev(f(:,3),3) = fcol;
                y            = max(ev')';
                
                % vertex interpolated colour
                mesh.h.FaceVertexCData = y;
                mesh.h.FaceColor = 'interp';
                
                data.overlay.vertexcdata = y;
                data.overlay.facecdata   = fcol;
                
            case 'vertex'
                % Set vertex color, using interpolated face colours
                fcol  = spm_mesh_smooth(mesh, fcol(:), 4);
                fcol(isnan(fcol)) = 0;
                fcol  = S(1) + ((S(2)-S(1))).*(fcol - min(fcol))./(max(fcol) - min(fcol));
                set(mesh.h,'FaceVertexCData',fcol(:),'FaceColor','interp');
        end
        
        % Use symmetric colourbar and jet as defaults
        s = max(abs(fcol(:))); caxis([-s s]);
        colormap('jet');
        alpha 1;
        
        % Return the face colours
        data.overlay.data  = fcol(:);       % the functional vector
        data.overlay.steps = step;          % the depths at which searched
        data.overlay.hits  = hits;          % num hits / intersects at each depth
        data.overlay.cast  = UseFaceVertex; % whether computed for faces or vertices
        
        data.overlay.FaceNormals   = FaceNorm;
        data.overlay.FaceCentroids = FaceCent;
        data.overlay.FaceNormLines = FaceNormLine;
        
        waitbar(1,wb,'Complete');
        close(wb);
    
    
    case 'spheres' % this would be better called 'box' in its current form
        
        % This method places a box (boundaries) around a sphere inflated around each
        % vertex point (a 'trap window') by a fixed radius. Mesh points
        % within these bounds are assigned to this vertex
        %
        % The functional overlay vector returned in the overlay
        % substructure contains 1 value per VERTEX and faces colours are
        % interpolated
        %
        
        debugplot = 0;
        
        r  = (nv/length(pos))*1.3;      % radius - number of closest points on mesh
        r  = max(r,1);                  % catch when the overlay is over specified!
        OL = sparse(length(L),nv);      % this will be overlay matrix we average
        w  = linspace(.1,1,r);          % weights for closest points
        w  = fliplr(w);                 % 
        M  = zeros( length(x), nv);     % weights matrix: size(len(mesh),len(AAL))
        
        fprintf('Using inside-spheres search algorithm\n');
        tic
        for i = 1:length(x)
            if any(L(i))      
                newv = [];
                r   = 7;
                res = 20;
                th  = 0:pi/res:2*pi;
                r0  = [th(1:2:end-1) th(end) fliplr(th(1:2:end-1))];  
                
                % make [circle] radius change with z-direction (height)
                r0 = th.*fliplr(th);
                r0 = r0/max(r0);
                r0 = r0*r;
                r0 = r0 ;
                
                % the height at which each circle making the sphere will go
                z0 = linspace(v(i,3)-r,v(i,3)+r,(res*2)+1);

                % this generates the vertices of the sphere
                for zi = 1:length(z0)
                    xunit = r0(zi) * cos(th) + v(i,1);
                    yunit = r0(zi) * sin(th) + v(i,2);
                    zunit = repmat(z0(zi),[1,length(xunit)]);
                    newv  = [newv; [xunit' yunit' zunit']];
                end

                if debugplot
                    hold on;
                    s1 = scatter3(v(i,1),v(i,2),v(i,3),200,'r','filled');
                    s2 = scatter3(newv(:,1),newv(:,2),newv(:,3),150,'b');
                    s2.MarkerEdgeAlpha = 0.1;
                    drawnow;
                end

                % Determine whether this point if left or right hemisphere
                LR     = v(i,1);
                IsLeft = (LR-cnt(1)) < 0;
                
                if IsLeft; lri = Lft;
                else;      lri = Rht;
                end
                
                % Bounding box
                bx = [min(newv); max(newv)];
                inside = ...
                    [mv(lri,1) > bx(1,1) & mv(lri,1) < bx(2,1) &...
                     mv(lri,2) > bx(1,2) & mv(lri,2) < bx(2,2) &...
                     mv(lri,3) > bx(1,3) & mv(lri,3) < bx(2,3) ];
                 
                ind = lri(find(inside));
                OL(i,ind) = L(i);
                M (i,ind) = 1;
                indz{i}   = ind;
                w         = 1;
            end
        end
        stime = toc;
        fprintf('Routine took %d seconds\n',stime);

    case 'euclidean'
        
        % Computes (vectorised) euclidean distance from each vertex to
        % every mesh point. Selects closest n to represent vertex values and 
        % weights by distabnce. n is defined by nmeshpoints / nvertex *1.3
        %
        % The functional overlay vector returned in the overlay
        % substructure contains 1 value per VERTEX and face colours are
        % interpolated
        %
        
        debugplot = 0;
        
        OL = sparse(length(L),nv);      % this will be overlay matrix we average
        r  = (nv/length(pos))*1.3;      % radius - number of closest points on mesh
        r  = max(r,1);                  % catch when the overlay is over specified!
        w  = linspace(.1,1,r);          % weights for closest points
        w  = fliplr(w);                 % 
        M  = zeros( length(x), nv);     % weights matrix: size(len(mesh),len(AAL))

        fprintf('Using euclidean search algorithm\n');
        tic
        for i = 1:length(x)
            
            % Print progress
            if i > 1; fprintf(repmat('\b',[size(str)])); end
            str = sprintf('%d/%d',i,(length(x)));
            fprintf(str);

            % Restrict search to this hemisphere
            LR     = v(i,1);
            IsLeft = (LR-cnt(1)) < 0;
            
            if IsLeft; lri = Lft;
            else;      lri = Rht;
            end

            % Compute euclidean distances
            dist       = cdist(mv(lri,:),v(i,:));
            [junk,ind] = maxpoints(dist,max(r,1),'min');
            ind        = lri(ind);
            OL(i,ind)  = w*L(i);
            M (i,ind)  = w;
            indz(i,:)  = ind;

            if debugplot
                hold on;
                s1 = scatter3(v(i,1),v(i,2),v(i,3),200,'r','filled');
                s2 = scatter3(mv(ind,1),mv(ind,2),mv(ind,3),150,'b','filled');
                drawnow;
            end

        end
        stime = toc;
        fprintf('Routine took %d seconds\n',stime);
        
        
    case 'aal'
        
        % new method for AAL90: wrapper on ray casting routine
        load DenseAAL.mat

        ol    = zeros(length(v),1);
        for i = 1:length(L)
            these = find(vi==i);
            ol(these) = L(i);
        end
        
        % update sourcemodel
        data.sourcemodel.pos = v;
        data.overlay.orig    = ol;
        data.overlay.method  = data.overlay.method{2};
        
        data = aplot.overlay(data,ol,write,fname,colbar);
        return;
        
    case 'aal_light'
        
        % new method for AAL90: wrapper on ray casting routine
        load LightAAL.mat

        ol    = zeros(length(v),1);
        for i = 1:length(L)
            these = find(vi==i);
            ol(these) = L(i);
        end
        
        % update sourcemodel
        data.sourcemodel.pos = v;
        data.overlay.orig    = ol;
        data.overlay.method  = data.overlay.method{2};
        data.overlay.atlasvalues = L;
        data = aplot.overlay(data,ol,write,fname,colbar);
        return;
        
end


switch lower(method)
    case {'raycast','aal'}
        % Don't do anything        
    otherwise
        
        % kill interhems!
        VL = find(v(:,1) < 0);
        VR = find(v(:,1) > 0);
        ML = find(mv(:,1) < 0);
        MR = find(mv(:,1) > 0);
        OL(VL,MR) = 0;
        OL(VR,ML) = 0;
        
        
        fprintf('\n'); clear L;
        if ~interpl
             % mean value of a given vertex
            OL = mean((OL),1);
        else
            for i = 1:size(OL,2)
                % average overlapping voxels
                L(i) = sum( OL(:,i) ) / length(find(OL(:,i))) ;
                NumComp(i) =  length(find(OL(:,i)));
            end
            OL = L;
        end

        % normalise and rescale
        OL = double(full(OL));
        y  = S(1) + ((S(2)-S(1))).*(OL - min(OL))./(max(OL) - min(OL));

        y(isnan(y)) = 0;
        y  = full(y);
        y  = double(y);

        % spm mesh smoothing
        %--------------------------------------------------------------------------
        fprintf('Smoothing overlay...\n');
        y  = spm_mesh_smooth(mesh, y(:), 4);
        y(isnan(y)) = 0;
        y  = S(1) + ((S(2)-S(1))).*(OL - min(OL))./(max(OL) - min(OL));
        y(isnan(y)) = 0;

        % return these in data structre
        data.overlay.data           = y;
        data.overlay.smooth_weights = M;
        data.overlay.NumComp        = NumComp;
        data.overlay.indz           = indz;
        data.overlay.w              = w;

        set(mesh.h,'FaceVertexCData',y(:),'FaceColor','interp');
        drawnow;
        shading interp
        % force symmetric caxis bounds
        s = max(abs(y(:))); caxis([-s s]);
        colormap('jet');
        alpha 1;
end


if colbar
    data.overlay.cb = InteractiveColorbar;
end
    
% switches for writing out other file formats
%-------------------------------------------------------------------------
if write == 1;
    fprintf('Writing overlay gifti file: %s\n',[fname 'Overlay.gii']);
    g       = gifti;
    g.cdata = double(y);
    g.private.metadata(1).name  = 'SurfaceID';
    g.private.metadata(1).value = [fname 'Overlay.gii'];
    save(g, [fname  'Overlay.gii']);
elseif write == 2 
        % write STL
        fprintf('Writing mesh and overlay as STL object\n');
        m.vertices = double(mesh.vertices);
        m.faces    = double(mesh.faces);
        y = spm_mesh_smooth(mesh, y(:), 8); % hard smoothing
        y          = double(y);
        cdata      = mean(y(m.faces),2);
        
        I = cdata;
        Colors   = jet;
        NoColors = length(Colors);

        Ireduced = (I-min(I))/(max(I)-min(I))*(NoColors-1)+1;
        RGB      = interp1(1:NoColors,Colors,Ireduced);
        
        fCols8bit= RGB*255;
        stlwrite([fname '.stl'],m,'FaceColor',fCols8bit)
elseif write == 3 
    % write vrml
    fprintf('Writing vrml (.wrl) 3D object\n');
    vrml(gcf,[fname]);
elseif write == 4
    fprintf('Generating nifti volume for writing\n');    
    new = mesh;
    dim = ceil(nthroot(length(y),3));
    V   = sm2vol(new.vertices,dim*3,y,256);
    % just return the volume in the output for now
    data.overlay.volume = V;
end

end
drawnow;


% if a timefreqanal structure was passed: allow click-for-plot (see Mylcmv)
%--------------------------------------------------------------------------
if isfield(data.overlay,'tf_interactive') && isstruct(data.overlay.tf_interactive)
    tf   = data.overlay.tf_interactive;
    f0   = get(gca,'parent');
    f    = figure('position',[1531         560         560         420]);
    %waitfor(f) % while the box is open
    fprintf('Waiting: Click in table to view corresponding time-freq!\n');

    %while isvalid(f)
    H    = datacursormode(f0);
    H.DisplayStyle = 'datatip';

    ch = get(f0,'children');
    waitfor(f0.WindowButtonMotionFcn)
            H    = datacursormode(f0);
            INFO = getCursorInfo(H);
            
            if isstruct(INFO) && isfield(INFO,'Position')
                current      = INFO.Position;
                fprintf('Got point: fetching time-freq plot\n');
                [this,thisi] = find( cdist(current,INFO.Target.Vertices) == 0);
                
                figure(f);
                tf.index = thisi; % choose vertex (source) index to plot, rerun for plot
                timefreqanal(tf);
                drawnow;
            end
    %end
end

% if a spatial PCA was requested
%--------------------------------------------------------------------------
if isfield(data.overlay,'pca')
    if data.overlay.pca
        f  = mesh.faces;
        A  = spm_mesh_adjacency(f);
        sy = double(y(:));
        sy = sy - mean(abs(sy(:))); 
        ya = sy.*A;
        
        % hardwire an assumption that you wouldnt look for more than 4
        pks = findpeaks(y);
        pks = min(length(pks),4); 
        
        %[U,S,V] = svds(ya,pks);
        %for pc = 1:length(S)
        %    comp(pc,:) = sparse(U(:,pc)*S(pc,pc)*mean(V(:,pc))');
        %end
        
        [comp,D] = eigs(ya,pks); % use eigs instead
        pc   = size(comp,2);
        comp = comp';
        comp = D*comp;
        
        ncomp = pc;
        fprintf('Found %d principal components\n',pc);
        data.overlay.pca = comp;
        
        f0 = get(gca,'parent');
        f  = figure('position',[1531         560         560         420]);
        t  = uitable(f);
        for i = 1:ncomp
            d{i,1} = sprintf('Component %d',i);
            d{i,2} = false;
        end
        d = [ {'All' , false }; d];

        t.Data = d;
        t.ColumnName = {'PC'};
        t.ColumnEditable = true;

        %waitfor(f) % while the peaks box is open
        fprintf('Waiting: Click in table to view components!\n');
        while isvalid(f)
            try
                waitfor(t,'Data')
                i = find(cell2mat(t.Data(:,2)));
                if any(i)
                    if i > 1
                        i  = i - 1;
                        if size(comp,2) == length(mesh.vertices)
                            this = full(comp(i,:));
                        else
                            this = full(comp(i,:)*M');
                        end
                        if all(size(this) > 1)
                            this = sum(this,1);
                        end

                        Y = spm_mesh_smooth(mesh, this(:), 8);

                        thefig = get(f0,'children');
                        hh = get(thefig(end),'children');
                        set(hh(end),'FaceVertexCData',Y(:),'FaceColor','interp');
                        set(hh(end),'FaceAlpha',1);
                        drawnow;
                        shading interp
                    else
                        %otherwise just plot the whole lot
                        thefig = get(f0,'children');
                        hh = get(thefig(end),'children');
                        set(hh(end),'FaceVertexCData',y(:),'FaceColor','interp');
                        drawnow;
                        shading interp
                    end
                end
            catch
                return;
            end
        
        end  
    end
end

% if search of local maxima was requested
%--------------------------------------------------------------------------
if isfield(data.overlay,'components')
    if data.overlay.components
                
        fprintf('Computing local maxima on surface\n');
        if length(data.overlay.orig) ~= length(mesh.vertices);
            msh.vertices = data.sourcemodel.pos;
            msh.faces    = delaunayTriangulation(msh.vertices);
            msh.faces    = msh.faces(:,1:3);
            useM = 1;
        else
            msh  = mesh;
            useM = 0;
        end
        A   = spm_mesh_adjacency(msh);
        T   = data.overlay.orig;
        out = isnan(T(:)');
        Lm  = [];
        for i=find(~out)
            v = T(logical(A(i,:)));
            if ~any(v>T(i))
                Lm = [Lm i];
            end
        end
        
        % sort these by size
        fprintf('Sorting...\n');
        maxs   = T(Lm);
        [~,Ii] = sort(abs(maxs),'descend');
        P      = Lm(Ii);
        ncomp  = length(P);
        
        
        f0 = get(gca,'parent');
        f  = figure('position',[1531         560         560         420]);
        t  = uitable(f);
        for i = 1:ncomp
            d{i,1} = sprintf('Comp %d',i);
            d{i,2} = false;
        end
        d = [ {'All' , false }; d];

        t.Data = d;
        t.ColumnName = {'Component'};
        t.ColumnEditable = true;

        %waitfor(f) % while the peaks box is open
        fprintf('Waiting: Click in peaks table to view peaks!\n');
        while isvalid(f)

            waitfor(t,'Data')
            i = find(cell2mat(t.Data(:,2)));
            if any(i)
                if i > 1
                    i  = i - 1;
                    if useM
                        dM = M*0;
                        dM(P(i),:) = M(P(i),:);
                    else
                        dM = zeros(length(T),1);
                        dM(P(i)) = 1;
                    end
                    try    this = data.overlay.orig*dM;
                    catch; this = dM*T;
                    end
                    Y = spm_mesh_smooth(mesh, this(:), 4);

                    thefig = get(f0,'children');
                    hh = get(thefig(end),'children');
                    set(hh(end),'FaceVertexCData',Y(:),'FaceColor','interp');
                    drawnow;
                    shading interp
                else
                    %otherwise just plot the whole lot
                    thefig = get(f0,'children');
                    hh = get(thefig(end),'children');
                    set(hh(end),'FaceVertexCData',y(:),'FaceColor','interp');
                    drawnow;
                    shading interp
                end
            end

        end  
    end
end

% If 'Peaks' was requested while using an AAL atlas
%--------------------------------------------------------------------------
if isfield(data.overlay,'Peaks')
    if isfield(data.overlay.Peaks,'Labels')
        f0 = get(gca,'parent');
        f = figure('position',[1531         560         560         420]);
        t = uitable(f);
        for i = 1:length(data.overlay.Peaks.Labels)
            d{i,1} = data.overlay.Peaks.Labels{i};
            d{i,2} = data.overlay.Peaks.Values(i);
            d{i,3} = false;
        end
        d{i+1,1} = 'All';
        d{i+1,2} = '--';
        d{i+1,3} = false;
        
        t.Data = d;
        t.ColumnName = {'Position','Val','Spotlight'};
        t.ColumnEditable = true;
        
        %waitfor(f) % while the peaks box is open
        fprintf('Waiting: Click in peaks table to view peaks!\n');
        while isvalid(f)
            
            waitfor(t,'Data')
            i = find(cell2mat(t.Data(:,3)));
            if any(i)
                if i < length(d)
                    this = t.Data(i,:);
                    % work backwards to project only this component
                    thislab = find(strcmp(this{1},data.atlas.AAL_Labels));
                    dM      = M;
                    n       = 1:size(dM,1);
                    dM(find(~ismember(n,thislab)),:) = 0;
                    dM = dM'*data.overlay.orig(:);
                    dM = full(double(dM));
                    Y = spm_mesh_smooth(mesh, dM(:), 4);
                    
                    thefig = get(f0,'children');
                    hh = get(thefig(end),'children');
                    set(hh(end),'FaceVertexCData',Y(:),'FaceColor','interp');
                    set(hh(end),'FaceAlpha',1);
                    drawnow;
                    shading interp
                else
                    % otherwise just plot the whole lot
                    thefig = get(f0,'children');
                    hh = get(thefig(end),'children');
                    set(hh(end),'FaceVertexCData',y(:),'FaceColor','interp');
                    drawnow;
                    shading interp
                end
            end
        
        end
        
    end
end

end
