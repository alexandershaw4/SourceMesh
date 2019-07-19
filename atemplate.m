function data = atemplate(varargin)
% Add networks and overlays to a smoothed brain mesh or gifti object.
%
% Plots: volumes, surfaces, overlays, networks, labels, videos. 
% Examples and usages below.
%
% If you get errors using the mex files, delete them. 
% 
% Now uses fieldtrip to read nifti's.
%
% MESHES:
%--------------------------------------------------------------------------
%
%  % Plot the default template mesh:
%  atemplate()         
%
%  % Plot a supplied (gifti) mesh:
%  atemplate('gifti',mesh)   
%
%  % Plot mesh & write out gifti:
%  atemplate('gifti',mesh, 'write',name);  
%  
%  % Plot mesh from nifti volume:
%  atemplate('mesh','mymri.nii')
%
%  % Plot only one hemisphere:
%  atemplate('hemi','left'); atemplate('hemi','L'); atemplate('hemi','l');
%  atemplate('hemi','right');atemplate('hemi','R'); atemplate('hemi','r');
%  atemplate('gifti',mesh,'hemi','left') 
%
%  % Fill holes that appear from hemisphere separation
%  atemplate('gifti',mesh,'hemi','left','fillholes') 
%
%  Plot volume and supply affine transformation matrix:
%  atemplate('mesh','mymri.nii','affine',affinematrix)
%
%
% OVERLAYS:
%--------------------------------------------------------------------------
%
%  % Plot template mesh with overlay from AAL90. L is [90x1]
%  atemplate('overlay',L,'method,{'aal','spheres'});   
%
%  % Plot overlay aligned to mesh using a Euclidean ICP search:
%  atemplate('overlay',L,'sourcemodel',pos,'method','euclidean');  
%
%  % Plot overlay aligned to mesh using a sphere-based trap radius method (def):
%  atemplate('overlay',L,'sourcemodel',pos,'method','spheres');  
%
%  % Plot overlay aligned to mesh using ray casting
%  atemplate('overlay',L,'sourcemodel',pos,'method','raycast');  
%  atemplate('overlay',L,'sourcemodel',pos,'method','raycast','depth',-1.5:.05:1.5);  
%
%  % Plot template with overlay values L at sourcemodel values sormod, interpolated on surface.
%  % Sormod is n-by-3, L is n-by-1.
%  atemplate('sourcemodel',sormod,'overlay',L)  
%
%  % Plot the supplied gifti mesh with overlay values L at sourcemodel locations 
%  % sormod interpolated on surface. 
%  % Sormod is n-by-3, L is n-by-1.
%  atemplate('gifti',mesh,'sourcemodel',sormod,'overlay',L)  
%
%  %  - Plot as above but write out TWO gifti files:
%  %  1. MYGifti.gii is the gifti mesh 
%  %  2. MYGiftiOverlay.gii is the corresponding overlay data
%  atemplate('gifti',mesh,'sourcemodel',sormod,'overlay',L,'write','MYGifti')  
%
%  % Return nifti volume:
%  atemplate('gifti',mesh,'sourcemodel',sormod,'overlay',L,'writenii','mynifti')
%
%  % Plot overlay from nifti volume
%  atemplate('overlay','overlay_volume.nii')
%
%
% % Plot overlay (e.g. from nifti volume) along with curvature on an
% inflated mesh:
%  afigure;D=atemplate('mesh','def2','inflate',...
%                      'overlay',{'curvature','Vis100B_Av.nii'},...
%                      'thresh',.4,'open');
%
%
%  *Note on sourcemodel option: Some fieldtrip sourcemodels have x & y
%  swapped (?), undo by doing sm = [sm(:,2),sm(:,1),sm(:,3)];
%
%  % Co-register the surfaces of the nii volumes in mesh and overlay,
%  % put in aal90 space and add labels:
%  atemplate('mesh',t1.nii,'overlay',functional.nii,'template','aal90','labels')
%
%
%  % Put overlay in AAL space and use interactive 'peaks' (clickable)
%  atemplate('sourcemodel',sormod,'overlay',randi([0 9],1000,1),'template','aal90','peaks')
%
%
%  % Find local maxima in overlay:
%  atemplate('sourcemodel',sormod,'overlay',randi([0 9],1000,1),'components','nocolbar')
%
%  % Do PCA on overlay
%  atemplate('sourcemodel',sormod,'overlay',randi([0 9],1000,1),'pca','nocolbar')
%
%
%
% VIDEO OVERLAY:
%--------------------------------------------------------------------------
%
%  % Plot a video overlay and write it out:
%  atemplate('gifti',g,'sourcemodel',sormod,'video',m,'name',times); 
%
%  % Where:
%  - g      = the gifti surface to plot
%  - sormod = sourcemodel vertices
%  - m      = overlay values [vertices * ntimes] 
%  - name   = video savename
%  - times  = vector of titles (time values?)
%
%
% NETWORKS:
%--------------------------------------------------------------------------
%
%  % Plot template mesh with 90x90 AAL network, A:
%  atemplate('network',A); 
%
%  % Plot network A  at sourcemodel locations in 'sormod'. 
%  % Sormod is n-by-3, network is n-by-n.
%  atemplate('sourcemodel',sormod,'network',A);  
%
%  % As above but writes out .node and .edge files for the network, and the gifti mesh file.
%  atemplate('sourcemodel',sormod,'network',A,'write','savename'); 
%  
%  % Plot network defined by .edge and .node files*:
%  atemplate('network','edgefile.edge');
%  % Note, this option sets the 'sourcemodel' coordinates to the vertices
%  defined in the .node file, unless flag to register to atlas space
%
%
% Project to ATLAS
%--------------------------------------------------------------------------
%
%  % Put overlay into atlas space: [choose aal90, aal78 or aal58]
%  atemplate('sourcemodel',sormod,'overlay',o,'template','aal58')
%
%  % Put network into atlas space: 
%  atemplate('sourcemodel',sormod,'network',N,'template','aal78')
%
%  % Put video into atlas space: 
%  atemplate('sourcemodel',sormod,'video',m,'name',times,'template','aal78')
%
%
% OTHER:
%--------------------------------------------------------------------------
%
%  % Export 3D images (overlays, meshes, networks) as VRML & .stl:
%  atemplate( ... ,'writestl','filename.stl');
%  atemplate( ... ,'writevrml','filename.wrl');
%
%
%  % Plot default AAL90 node labels on default mesh:
%  atemplate('labels');         
%
%  % Plot specified labels at centre of roi's specified by all_roi_tissueindex:
%  atemplate('labels', all_roi_tissueindex, labels); 
%
%  % Where:
%  % all_roi_tissue = a 1-by-num-vertices vector containing indices of the
% roi this vertex belongs to
%  % 'labels' = the labels for each roi. 
%  % The text labels are added at the centre of the ROI.
%  
%  Labels notes:
%     - If plotting a network, only edge-connected nodes are labelled.
%     - If plotting a set of nodes (below), only those are labelled.
%     - Otherwise, all ROIs/node labels are added!
%
%  % Plot dots at node==1, i.e. N=[90,1]:
%  atemplate('nodes', N);             
%
%  % Plot tracks loaded with trk_read, from along-tract-stats toolbox.
%  % This function requires some work...
%  atemplate('tracks',tracks,header); 
%
%  Any combination of the inputs should be possible.
%  See scripts in 'Examples' folder for more help.
%
% 
%
%
%
% AN EXAMPLE NETWORK [1]: from 5061 vertex sourcemodel with AAL90 labels
%--------------------------------------------------------------------------
% load New_AALROI_6mm.mat          % load ft source model, labels and roi_inds
% net  = randi([0 1],5061,5061);   % generate a network for this sourmod
% pos  = template_sourcemodel.pos; % get sourcemodel vertices
% labs = AAL_Labels;               % roi labels
% rois = all_roi_tissueindex;      % roi vertex indices
%
% atemplate('sourcemodel',pos,'network',net,'labels',rois,labs);
%
%
% AN EXAMPLE NETWORK [2]: from volume and node/edge files, put in aal58 space:
%--------------------------------------------------------------------------
% atemplate('mesh',t1.nii,'network','test_sourcemod.edge','template','aal58')
%
%
% Requirements: SPM (for spm_vec.m), fieldtrip (for ft_read_mri.m)
%
% AS17


% defaults
%--------------------------------------------------------------------------
data         = struct;     % main output structure
in.pmesh     = 1;          % flag to display the mesh or not
in.labels    = 0;          % flag to display labels on rois
in.write     = 0;          % flag to write mesh & overlay
in.fname     = [];         % filename for saves (e.g. for sae gifti)
in.fighnd    = [];         % put the plot into an existing figure handle
in.colbar    = 1;          % display colourbar
in.template  = 0;          % register source model to a template using ICP
in.orthog    = 0;          % orthogonalise [ignore]
in.inflate   = 0;          % inflate the mesh [just calls spm_mesh_inflate]
in.peaks     = 0;          % compute peaks on surface functional overlay
in.components = 0;         % compute local maxima on functional overlay
in.pca        = 0;         % compute pca on surface functional overlay
in.flip       = 0;         % ignore
in.affine     = 0;         % supply affine transformation matrix
in.netcmap    = 0;         % colormap for network
in.method     = 'raycast'; % volume to surface algorithm
in.depth      = [];        % depth for raycasting
in.all_roi_tissueindex = []; % indices of roi each voxel belongs to
in.thelabels  = [];        % supply cell array of labels for rois
in.tf_interactive = 0;     % ignore
in.checkori       = 0;     % interactively check/adjust orientation of mesh
in.fillholes      = 0;     % fill holes in mesh 
in.hemi = 'both';          % hemisphere to plot
in.roi  = [];
in.optimise = 0;
in.post_parcel = [];
in.thresh      = [];
in.open        = 0;

% specified inputs [override defaults]
%--------------------------------------------------------------------------
for i  = 1:length(varargin)
    if strcmp(varargin{i},'overlay');     in.L   = varargin{i+1}; end
    if strcmp(varargin{i},'open');        in.open= 1;             end
    if strcmp(varargin{i},'hemi');        in.hemi= varargin{i+1}; end
    if strcmp(varargin{i},'peaks');       in.peaks = 1;           end
    if strcmp(varargin{i},'sourcemodel'); in.pos = varargin{i+1}; end
    if strcmp(varargin{i},'roi');         in.roi = varargin{i+1}; end
    if strcmp(varargin{i},'network');     in.A   = varargin{i+1}; end
    if strcmp(varargin{i},'tracks');      in.T   = varargin{i+1}; in.H = varargin{i+2}; end
    if strcmp(varargin{i},'nosurf');      in.pmesh  = 0;            end
    if strcmp(varargin{i},'nodes');       in.N = varargin{i+1};     end
    if strcmp(varargin{i},'gifti');       in.g = varargin{i+1};     end
    if strcmp(varargin{i},'mesh');        in.g = varargin{i+1};     end
    if strcmp(varargin{i},'optimise');    in.optimise = varargin{i+1};end
    if strcmp(varargin{i},'post_parcel'); in.post_parcel = varargin{i+1};end
    if strcmp(varargin{i},'fillholes');   in.fillholes = 1;         end
    if strcmp(varargin{i},'thresh');      in.thresh = varargin{i+1};end
    if strcmp(varargin{i},'affine');      in.affine = varargin{i+1};end
    if strcmp(varargin{i},'funcaffine');  in.funcaffine = varargin{i+1}; end
    if strcmp(varargin{i},'inflate');     in.inflate = 1;           end
    if strcmp(varargin{i},'orthog');      in.orthog = varargin{i+1};end
    if strcmp(varargin{i},'components');  in.components = 1;        end
    if strcmp(varargin{i},'pca');         in.pca = 1;               end    
    if strcmp(varargin{i},'flip');        in.flip = 1;              end
    if strcmp(varargin{i},'checkori');    in.checkori = 1;          end
    if strcmp(varargin{i},'netcmap');     in.netcmap= varargin{i+1};end
    if strcmp(varargin{i},'write');       in.write  = 1; in.fname = varargin{i+1}; end
    if strcmp(varargin{i},'writestl');    in.write  = 2; in.fname = varargin{i+1}; end
    if strcmp(varargin{i},'writevrml');   in.write  = 3; in.fname = varargin{i+1}; end
    if strcmp(varargin{i},'writenii');    in.write  = 4; in.fname = varargin{i+1}; end
    if strcmp(varargin{i},'fighnd');      in.fighnd = varargin{i+1}; end
    if strcmp(varargin{i},'nocolbar');    in.colbar = 0;             end
    if strcmp(varargin{i},'method');      in.method = varargin{i+1}; end
    if strcmp(varargin{i},'depth');       in.depth  = varargin{i+1}; end
    if strcmp(varargin{i},'video');       in.V     = varargin{i+1}; 
                                          in.fpath = varargin{i+2}; 
                                          in.times = varargin{i+3}; end
    if strcmp(varargin{i},'othermesh');   in.M = varargin{i+1}; in.O = varargin{i+2};   end  
    if strcmp(varargin{i},'tf_interactive');in.tf_interactive = varargin{i+1}; end
    if strcmp(varargin{i},'labels');      in.labels = 1;
        try in.all_roi_tissueindex = varargin{i+1};
            in.thelabels = varargin{i+2};
        end
    end
    if strcmp(varargin{i},'template')
        in.template = 1;
        in.model    = varargin{i+1};
    end  
    
    % Allow passing of existing atemplate-returned structure 
    if isstruct(varargin{i}) && isfield(varargin{i},'in')
        fprintf('User specified plot structure\n');
        data = varargin{i};
        mesh = parse_mesh(data.mesh,data.in,data);
        data = parse_plots(data,data.in);
        return
    end
end



% Preliminaries & parse functions for plots
%--------------------------------------------------------------------------
if in.open
    data = isopen(data,in);
    return;
end

data = sort_sourcemodel(data,in);       % Sourcemodel vertices

[mesh,data] = get_mesh(in,data);        % Get Surface

[data,in] = sort_template(data,in);     % Template space? (aal90/78/58)

[mesh,data] = parse_mesh(mesh,in,data); % Mesh to put stuff on
data.mesh   = mesh;

data = parse_plots(data,in); % Overlays, networks, tracts, nodes, videos etc.

data.in = in; % Return the input options for re-run

end




% FUNCTIONS
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function data = isopen(data,in)
    % a wrapper on atemplate for generating the same plot for the left and
    % right hemis separately in different sublots
    s(1) = subplot(121);
    s(2) = subplot(122);
    
    data_in     = data;
    
    % Plot 1.
    in.hemi     = 'l';
    in.fighnd   = s(1);
    data        = sort_sourcemodel(data_in,in); % Sourcemodel vertices
    [mesh,data] = get_mesh(in,data);         % Get Surface
    [data,in]   = sort_template(data,in);    % Template space? (aal90/78/58)
    [mesh,data] = parse_mesh(mesh,in,data);  % Mesh to put stuff on
    data.mesh   = mesh;
    data        = parse_plots(data,in);      % Overlays, networks, etc
    data.in     = in;                        % Return everything
    data_l      = data;
    
    view([90 0]);
    
    % Plot 2.
    in.hemi     = 'r';
    in.fighnd   = s(2);
    data        = sort_sourcemodel(data_in,in); % Sourcemodel vertices
    [mesh,data] = get_mesh(in,data);         % Get Surface
    [data,in]   = sort_template(data,in);    % Template space? (aal90/78/58)
    [mesh,data] = parse_mesh(mesh,in,data);  % Mesh to put stuff on
    data.mesh   = mesh;
    data        = parse_plots(data,in);      % Overlays, networks, etc
    data.in     = in;                        % Return everything
    data_r      = data;
    
    view([270 0]);
    
    data   = [];
    data.l = data_l;
    data.r = data_r;

end

function data = parse_plots(data,i)

% unpack triggers
inputs = i;

% overlays
if isfield(inputs,'L')
    % copy over overlay options
    data.overlay.peaks      = i.peaks;
    data.overlay.components = i.components;
    data.overlay.pca        = i.pca;
    data.overlay.method     = i.method;
    data.overlay.depth      = i.depth;
    data.overlay.thresh     = i.thresh;
    data.overlay.tf_interactive = i.tf_interactive;
    
    if isfield(i,'funcaffine')
        data.overlay.affine = i.funcaffine;
    end
    
    data = overlay(data, (i.L),i.write,i.fname,i.colbar);
end 

isover = exist('L','var') || exist('V','var');
if  isover && exist('A','var') 
    i.colbar = 0;
    alpha(.2);
end

% networks
if isfield(inputs,'A')
    data = connections(data,i.A,i.colbar,i.write,i.fname,i.netcmap); 
end 

% tracts
if isfield(inputs,'T')
    data = drawtracks(data,i.T,i.H);                  
end 

% nodes
if isfield(inputs,'N')
    data = drawnodes(data, i.N);                 
end 

% labels
data = parse_labels(i,data);

% video
if isfield(inputs,'V')
    tv = 1:size(i.V,2);
    try tv = i.times; end
    data = video(data,i.V,1,i.fpath,tv); 
end


end

function data = parse_labels(i,data)
% decide which labels to include depending on what we're plotting
if i.labels 
    if     isfield(i,'A')
                if isnumeric(i.A)
                    data = addlabels(data,i.A,i.all_roi_tissueindex,i.thelabels);
                elseif ischar(i.A)
                    E = data.network.edge;
                    data = addlabels(data,E,i.all_roi_tissueindex,i.thelabels);
                end
           
    elseif isfield(i,'N')
        if sum(ismember(size(i.N),[1 90])) == 2
            data = addlabels(data, diag(i.N),i.all_roi_tissueindex,i.thelabels);
        elseif sum(ismember(size(i.N),[1 90])) == 1
            data = addlabels(data, diag(sum(i.N,2)),i.all_roi_tissueindex,i.thelabels);
        end
        
    else;  n    = length(data.sourcemodel.pos);
           data = addlabels(data, ones(n,n),i.all_roi_tissueindex,i.thelabels);
    end
end
end

function data = sort_sourcemodel(data,i)
% Sort out what source model vertices we're going to use

if      isfield(i,'pos')
            fprintf('Using supplied sourcemodel vertices\n');
            pos = i.pos;
        
elseif  isfield(i,'A') && ischar(i.A)
        fprintf('Using coords in node-file as sourcemodel\n');
        
        [~,pos] = rw_edgenode(i.A); 
        pos = pos(:,1:3);
        
%  elseif  isfield(i,'L') && ischar(i.L)
%          %IF USING A NIFTI OVERLAY AS SOURCEMODEL, NEED THIS *before* TRYING TO DO
%          %TEMPLATE SPACE CONVERSION!
%          
%          [mesh,data] = get_mesh(i,data);
%          data.mesh = mesh;
%          [~,data]  = parse_overlay(i.L,data);
%          pos       = data.sourcemodel.pos;
         
else
        fprintf('Assuming AAL90 source vertices by default\n');
        load('AAL_SOURCEMOD');
        pos  = template_sourcemodel.pos;
end

if iscell(pos)
    % store 2nd input - vector descrbinig which vertices belong to which rois
    data.sourcemodel.vi = pos{2};
    pos                 = pos{1};
end

% Centre sourcemodel
dpos = pos - repmat(spherefit(pos),[size(pos,1),1]);

% Ensure stability if approaching centre (very tiny values)
if ~any(isnan(dpos(:)))
    pos = dpos;
end

data.sourcemodel.pos = pos;


% not hugely revelant here, but use this func to also store data for the
% post-hoc atlas-registration function if required
if iscell(i.post_parcel)
    data.post_parcel = i.post_parcel;
end


end

function [mesh,data] = parse_mesh(mesh,i,data)
% Figure out whether we actually want to plot a glass brain mesh, or not

% if only one hemisphere plot
hemi      = i.hemi;
data.hemi = hemi;

% if affine supplied or flip flag
affine = i.affine;
flip   = i.flip;

% if inflate, pass flag
inflate = i.inflate;

% check orientation?
checkori = i.checkori;

% fill holes during hemisphere separation?
dofillholes = i.fillholes;

% explicitly optimise the alignment of the cloud points?
optimise = i.optimise;

if     i.pmesh && ~isfield(i,'T')
       [mesh,data.sourcemodel.pos,h,p] = meshmesh(mesh,i.write,i.fname,i.fighnd,...
           .3,data.sourcemodel.pos,hemi,affine,flip,inflate,checkori,dofillholes,optimise);
elseif i.pmesh
       [mesh,data.sourcemodel.pos,h,p] = meshmesh(mesh,i.write,i.fname,i.fighnd,...
           .3,data.sourcemodel.pos,hemi,affine,flip,inflate,checkori,dofillholes,optimise);
else
    h = [];
    p = [];
end

mesh.h = h;
mesh.p = p;
end

function [mesh, data] = convert_mesh(mesh,data)
% Convert a STRUCTURAL volume to a surface mesh using an isosurface. Take
% some liberties with downsampling and rescaling 

if ischar(mesh)
    [fp,fn,fe] = fileparts(mesh);
    
    switch fe
        
        case{'.gz'}
            fprintf('Unpacking .gz\n');
            gunzip(mesh);
            mesh = strrep(mesh,'.gz','');
            [mesh, data] = convert_mesh(mesh,data);
            return;
            
        case{'.nii','.mri'}
            wb = waitbar(0,'Reading nifti volume...');
            
            % mini switch for load type
            switch fe
                case '.nii'
                    % load nifti volume file
                    fprintf('Reading Nifti volume\n');
                    
%                     ni    = load_nii(mesh);
%                     vol   = ni.img;
%                     
%                     % recover rotation matrix
%                     b = ni.hdr.hist.quatern_b;
%                     c = ni.hdr.hist.quatern_c;
%                     d = ni.hdr.hist.quatern_d;
%                     x = ni.hdr.hist.qoffset_x;
%                     y = ni.hdr.hist.qoffset_y;
%                     z = ni.hdr.hist.qoffset_z;
%                     a = sqrt(1.0-(b^2+c^2+d^2));
%                     
%                     R =  [ a*a+b*b-c*c-d*d   2*b*c-2*a*d       2*b*d+2*a*c     ;
%                            2*b*c+2*a*d       a*a+c*c-b*b-d*d   2*c*d-2*a*b     ;
%                            2*b*d-2*a*c       2*c*d+2*a*b       a*a+d*d-c*c-b*b ];
%                     t = [x; y; z];

                    % just use fieldtrip
                    ni  = ft_read_mri(mesh);
                    vol = ni.anatomy; 
                    
                    R = ni.transform(1:3,1:3);
                    t = ni.transform(1:3,end);
                    
                    % real-world coordinates = 
                    %[x; y; z] = R * [i; j; k] + t
                    
                    [nx,ny,nz] = size(vol);
                    
                    % compute edges and linearly interp grid coordinates
                    E1 = R*[1 ;1 ;1 ]+t;
                    E2 = R*[nx;ny;nz]+t;
                    
                    xarr = linspace(E1(1),E2(1),nx);
                    yarr = linspace(E1(2),E2(2),ny);
                    zarr = linspace(E1(3),E2(3),nz);
                                        
                    % NEW: subsample volume by 2. This speeds up subsequent patch
                    % reduction enormously with at relatively little cost to
                    % resolution
                    fprintf('Subsampling structural volume\n');
                    [xg,yg,zg] = meshgrid(yarr,xarr,zarr);
                    [NX, NY, NZ, vol] = reducevolume(xg,yg,zg,vol,2);
                    xarr = linspace(xarr(1),xarr(end),size(vol,1));
                    yarr = linspace(yarr(1),yarr(end),size(vol,2));
                    zarr = linspace(zarr(1),zarr(end),size(vol,3));

                    [X,Y,Z] = meshgrid(yarr,xarr,zarr);
                    
                    
                case '.mri'
                    % load ctf mri file
                    fprintf('Reading CTF_MRI4 volume\n');
                    ni  = ft_read_mri(mesh,'dataformat','ctf_mri4');
                    vol = ni.anatomy;
                    vol = vol - mean(vol(:));                    
            end
            
            
            if ndims(vol) ~= 3
                fprintf('Volume has wrong number of dimensions!\nUsing default mesh\n');
                mesh = read_nv;
                return
            end
            
            % bounds:
            waitbar(0.2,wb,'Reading nifti volume: Extracting isosurface');
            fprintf('Extracting ISO surface\n');
            B   = [min(data.sourcemodel.pos); max(data.sourcemodel.pos)];
            
            try    fv  = isosurface(X,Y,Z,vol,0.5); % nifti's 
            catch; fv  = isosurface(vol,0.5);       % .mri's / other volumes
            end
            
            % swap x y
            v  = fv.vertices;
            v  = [v(:,2) v(:,1) v(:,3)];
            fv.vertices = v;
             
            % reduce vertex density
            waitbar(0.6,wb,'Reading nifti volume: Reducing density');
            fprintf('Reducing patch density\n');
            nv  = length(fv.vertices);
            count  = 0;
                        
            while nv > 60000
               fv    = reducepatch(fv, 0.5);
               nv    = length(fv.vertices);
               count = count + 1;
            end
            
            %if count > 0
            %    fprintf('Smoothing surface\n');
            %    fv.vertices = sms(fv.vertices,fv.faces,1,2);
            %end

            % print
            waitbar(0.9,wb,'Reading nifti volume: Rescaling');
            fprintf('Patch reduction finished\n');
            fprintf('Rescaling mesh to sourcemodel\n');
            
            v = fv.vertices;
            for i = 1:3
                v(:,i) = rescale(v(:,i),B(:,i));
            end
            
            % return scaled mesh
            mesh            = [];
            mesh.nifti      = [fn fe];
            mesh.faces      = fv.faces;
            mesh.vertices   = v;
            mesh.vol        = vol;
            data.mesh       = mesh;
            
            close(wb);
            
        case{'.gii'}
            % load the gifti
            gi   = gifti(mesh);
            mesh          = [];
            mesh.faces    = gi.faces;
            mesh.vertices = gi.vertices;
            data.mesh     = mesh;
        
        otherwise
            
            % allow shorthands call for default meshes
            %-----------------------------------------
            if strcmp(lower(mesh),'def') || strcmp(lower(mesh),'def1')
                % ICBM152 - smoothed
                nmesh         = read_nv;
                mesh          = [];
                mesh.vertices = nmesh.vertices;
                mesh.faces    = nmesh.faces;
            end
            if strcmp(lower(mesh),'def2')
                % Ch2
                nmesh         = gifti('BrainMesh_Ch2.gii');
                mesh          = [];
                mesh.vertices = nmesh.vertices;
                mesh.faces    = nmesh.faces;
            end
            if strcmp(lower(mesh),'def3')
                % mni152_2009
                nmesh         = load('mni152_2009','faces','vertices');
                mesh          = [];
                mesh.vertices = nmesh.vertices;
                mesh.faces    = nmesh.faces;
            end
            if strcmp(lower(mesh),'def4')
                nmesh         = gifti('BrainMesh_ICBM152.gii');
                mesh          = [];
                mesh.vertices = nmesh.vertices;
                mesh.faces    = nmesh.faces;
            end
    end
    
elseif isnumeric(mesh) && ndims(mesh)==3
           
            % bounds:
            wb = waitbar(0,'Reading nifti volume...');
            waitbar(0.2,wb,'Reading nifti volume: Extracting isosurface');
            
            % NEW: subsample volume by 2. This speeds up subsequent patch
            % reduction enormously with at relatively little cost to
            % resolution
            fprintf('Subsampling structural volume\n');
            [NX, NY, NZ, mesh] = reducevolume(mesh,2);

            fprintf('Extracting ISO surface\n');
            B   = [min(data.sourcemodel.pos); max(data.sourcemodel.pos)];
            fv  = isosurface(mesh,0.5);            
            
            % swap x y
            v  = fv.vertices;
            v  = [v(:,2) v(:,1) v(:,3)];
            fv.vertices = v;
             
            % reduce vertex density
            fprintf('Reducing patch density\n');
            nv  = length(fv.vertices);
            count  = 0;
            
            waitbar(0.6,wb,'Reading nifti volume: Reducing density');
            while nv > 60000
                fv    = reducepatch(fv, 0.5);
                nv    = length(fv.vertices);
                count = count + 1;
            end

            % print
            fprintf('Patch reduction finished\n');
            fprintf('Rescaling mesh to sourcemodel\n');
            
            waitbar(0.9,wb,'Reading nifti volume: Rescaling');
            v = fv.vertices;
            for i = 1:3
                v(:,i) = rescale(v(:,i),B(:,i));
            end
            
            % return scaled mesh
            mesh            = [];
            mesh.faces      = fv.faces;
            mesh.vertices   = v;
            data.mesh       = mesh;    
            
            close(wb);        
end

end

function [data,i] = sort_template(data,i)
% If specified a template model, register data to it and return splined data as
% well as weights

if ~isfield(i,'pos')
    i.pos = data.sourcemodel.pos;
end
try
    data.template.model  = i.model;
    data.template.labels = i.labels;
end
if i.template
    atlas = dotemplate(i.model);
    rois  = get_roi_centres(atlas.template_sourcemodel.pos,atlas.all_roi_tissueindex);
    
    atlas.template_sourcemodel.pos = rois;
    atlas = rmfield(atlas,'all_roi_tissueindex');
    
    reg = interp_template(data.sourcemodel,rois);
    atlas.M    = reg.M;
    data.atlas = atlas;
    NM         = atlas.M;
    
    % rescale so not change amplitudes
    m  = max(NM(:));
    NM = NM/m; 
    
    data.atlas.M = NM;
    
    % generate & apply a ROI mask if requested
    %----------------------------------------------------
    if isfield(i,'roi')
        
        % should be cell array of labels
        if ischar(i.roi)
            i.roi = {i.roi};
        end
        
        % check rois exist
        labs     = atlas.AAL_Labels;
        for iroi = 1:length(i.roi)
            ROI_ind(iroi) = find(strcmp(i.roi{iroi},labs));
            
            if isempty( ROI_ind(iroi) )
                fprintf('Requested ROI (%s) not found\n');
            end
        end
        
        % apply mask!
        newM = NM*0;
        newM(:,ROI_ind) = NM(:,ROI_ind);
        NM   = newM;
    end
    
    % update sourcemodel and labels
    data.sourcemodel = atlas.template_sourcemodel;
    if i.labels; i.thelabels = atlas.AAL_Labels; end
    
    % overlay data
    if isfield(i,'L')
        if isnumeric(i.L) && ndims(i.L) ~= 3
            S  = [min(i.L(:)) max(i.L(:))];
            NL = i.L(:)'*NM;
            L  = S(1) + ((S(2)-S(1))).*(NL - min(NL))./(max(NL) - min(NL));
            L(isnan(L))=0;
            i.L = L;
        end
    end
    
    % network
    if isfield(i,'A')
        if isnumeric(i.A)
            S  = [min(i.A(:)) max(i.A(:))];
            NL = NM'*i.A*NM;
            A  = S(1) + ((S(2)-S(1))).*(NL - min(NL(:)))./(max(NL(:)) - min(NL(:)));
            A(isnan(A)) = 0;
            i.A = A;
        end
    end
    
    % video data
    if isfield(i,'V')
        S  = [min(i.V(:)) max(i.V(:))];
        for j = 1:size(i.V,2) % over time points
            NL(:,j) = i.V(:,j)'*NM;
        end
        V  = S(1) + ((S(2)-S(1))).*(NL - min(NL))./(max(NL) - min(NL));
        V(isnan(V))=0;
        if orthog
            % dont use this
            V = symm_orthog(V);
        end
        V(isnan(V))=0;
        i.V = V;
    end
        
end

end

function [mesh,data] = get_mesh(i,data)
% Decide what brain we're actually using, return it

try   mesh = i.g;
      fprintf('Using user provided mesh\n');
      
      if ischar(mesh) || isnumeric(mesh)
          [mesh,data] = convert_mesh(mesh,data);
      end
      
catch mesh = read_nv();
      fprintf('(Using template brain mesh)\n');
end

% if i.inflate
%     try fprintf('Trying to inflate mesh\n');
%         dmesh.vertices = mesh.vertices;
%         dmesh.faces    = mesh.faces;
%         dmesh = spm_mesh_inflate(dmesh,100);
%         mesh.vertices = dmesh.vertices;
%         mesh.faces    = dmesh.faces;
%     catch
%         fprintf('Couldnt find spm_mesh_inflate: is SPM installed?\n');
%     end
% end


end

function atlas = dotemplate(model)
% Put dense sourcemodel into an atlas space using ICP and linear
% interpolation
%
%
%

if ischar(model)
    switch model
        case lower({'aal','aal90'});   load New_AALROI_6mm.mat
        case lower('aal58');           load New_58cortical_AALROI_6mm
        case lower('aal78');           load New_AALROI_Cortical78_6mm
        otherwise
            fprintf('Model not found.\n');
            return;
    end
elseif iscell(model)
    template_sourcemodel.pos = model{1};
    all_roi_tissueindex      = model{2};
    AAL_Labels = [];
end

atlas.AAL_Labels = AAL_Labels;
atlas.all_roi_tissueindex = all_roi_tissueindex;
atlas.template_sourcemodel = template_sourcemodel;

end

function atlas = interp_template(atlas,pos)
% The Euclidean/ICP routine for atlas registration

if length(atlas.pos) == length(pos)
    fprintf('Overlay and atlas Vectors already match!\n');
    atlas.M = eye(length(pos));
    return;
end

fprintf('Scanning points:\n');
M = zeros( length(atlas.pos), length(pos) )';
r = floor( length(atlas.pos) / length(pos) );%1;
w = 1;

dist  = cdist(pos,atlas.pos)';    
%for i = 1:length(atlas.pos)
for i = 1:length(pos)    
    if i > 1; fprintf(repmat('\b',[size(str)])); end
    str = sprintf('%d/%d',i,(length(pos)));
    fprintf(str);

    [junk,ind] = maxpoints(dist(:,i),r,'min');
    M (i,ind)  = w;
end
fprintf('\n');
atlas.M = M';

end

function y = symm_orthog(x)
% Efficient symmetric orthogonalisation method for matrices, using:
%
% y = x * real(inv(x' * x)^(1/2));
%
% AS

fprintf('\nOrthogonalising\n');
y = [x] * real(inv([x]' * [x])^(1/2));
y = (max(x(:)) - min(x(:))) * ( (y - min(y(:))) / (max(y(:)) - min(y(:))) );

end

function data = connections(data,A,colbar,write,fname,netcmap)
% Network (Node & Edges) plotter function
%
%
pos = data.sourcemodel.pos;

if isfield(data.sourcemodel,'vi')
    % if pos is cell, it's because we've passed both a vertex set and a
    % vector that describes which vertex belongs to which roi
    v  = pos;
    vi = data.sourcemodel.vi;
    n  = unique(vi);
    roi = zeros(max(n),3);
    for i = 1:length(n)
        these = find(vi==n(i));
        roi(n(i),:) = spherefit(v(these,:));
    end
    pos = roi;
    for ip = 1:length(pos)
        if any(pos(ip,:))
            [~,this]  = min(cdist(pos(ip,:),data.mesh.vertices));
            pos(ip,:) = data.mesh.vertices(this,:);
        end
    end
end


% Read edge/node files if string
%--------------------------------------------------------------------------
if ischar(A)
    [fp,fn,fe]  = fileparts(A);
    [edge,node] = rw_edgenode(fn);
    A           = edge;
    
    if isfield(data,'template')
        if isfield(data.template,'model')
           fprintf('Doing atlas registration\n');
           i.template = 1;
           i.model    = data.template.model;
           i.labels   = data.template.labels;
           i.A        = A;
           data.sourcemodel.pos = node(:,1:3);
           [data,i]   = sort_template(data,i);
           A          = i.A;
        end
    end
end

A(isnan(A)) = 0;
A(isinf(A)) = 0;

% rescale network positions inside box boundaries of mesh
% (i thought meshmesh had already done this?)
if ~isempty(data.mesh.h)
    bounds = [min(data.mesh.vertices); max(data.mesh.vertices)];
    offset = 0.99;
    for ip = 1:3
        pos(:,ip) = bounds(1,ip) + ((bounds(2,ip)-bounds(1,ip))) .* ...
                    (pos(:,ip) - min(pos(:,ip)))./(max(pos(:,ip)) - min(pos(:,ip)));
        pos(:,ip) = pos(:,ip)*offset;
    end

    % redirect to closest mesh point (vertex)
    for ip = 1:length(pos)
        [~,this]  = min(cdist(pos(ip,:),data.mesh.vertices));
        pos(ip,:) = data.mesh.vertices(this,:);
    end
end

% Edges
%--------------------------------------------------------------------------
[node1,node2,strng] = matrix2nodes(A,pos);

% place both signed absmax value in overlay so that colorbar is symmetrical
strng2 = [strng; -max(abs(strng)); max(abs(strng))];
RGB    = makecolbar(strng2,netcmap);

% LineWidth (scaled) for strength
if any(strng)
    R = [-max(abs(strng)),max(abs(strng))];
    S = ( abs(strng) - R(1) ) + 1e-3;
    
    % If all edges same value, make thicker
    if  max(S(:)) == 1e-3; 
        S = 3*ones(size(S)); 
    end
else
    S = [0 0];
end

% If too few strengths, just use red edges
%--------------------------------------------------------------------------
LimC = 1;
if all(all(isnan(RGB)))
    RGB  = repmat([1 0 0],[size(RGB,1) 1]);
    LimC = 0;
end

data.network.edge = A;
data.network.node = pos;
data.network.RGB  = RGB;
data.network.tofrom.node1 = node1;
data.network.tofrom.node2 = node2;

% force girth boundaries if stable
if ~any(isnan( (S - min(S)) ./ (max(S) - min(S)) ))
    S = 0.1 + (3 - 0) .* (S - min(S)) ./ (max(S) - min(S));
end

% Paint edges
%--------------------------------------------------------------------------
for i = 1:size(node1,1)
    l0(i)=line([node1(i,1),node2(i,1)],...
        [node1(i,2),node2(i,2)],...
        [node1(i,3),node2(i,3)],...
        'LineWidth',S(i),'Color',[RGB(i,:)]);
end

% Set colorbar only if there are valid edges
%--------------------------------------------------------------------------
if any(i) && colbar
    set(gcf,'DefaultAxesColorOrder',RGB)
    set(gcf,'Colormap',RGB)
    if colbar
        %colormap(jet)
        %colorbar
        drawnow; pause(.5);
        a1  = gca;
        axb = axes('position', get(a1, 'position'));
        set(axb,'visible','off')
        axes(axb);
        %set(a1,'DefaultAxesColorOrder',RGB)
        set(gcf,'Colormap',RGB)
        
        if any(any(netcmap ~= 0)); 
                    colormap(netcmap);
        else;       colormap(jet);
        end
        
        colorbar('peer',a1,'South');
    end
end
if LimC && colbar
    axes(a1);
    caxis(R);
end

drawnow;


% Nodes (of edges only)
%--------------------------------------------------------------------------
% hold on;
% for i = 1:size(node1,1)
%     scatter3(node1(i,1),node1(i,2),node1(i,3),'filled','k');
%     scatter3(node2(i,1),node2(i,2),node2(i,3),'filled','k');
% end

drawnow;

if write;
   fprintf('Writing network: .edge & .node files\n');
   conmat2nodes(A,fname,'sourcemodel',pos);
end


end

function [node1,node2,strng] = matrix2nodes(A,pos)
% Write node & edge files for the AAL90 atlas
% Also returns node-to-node coordinates for the matrix specified.
%
% Input is the n-by-n connectivity matrix
% Input 2 is the sourcemodel vertices, n-by-3
%
% AS2017



node1 = []; node2 = []; strng = [];
for i = 1:length(A)
    [ix,iy,iv] = find(A(i,:));
    
    if ~isempty(ix)
        conns = max(length(ix),length(iy));
        for nc = 1:conns
            node1 = [node1; pos(i(1),:)];
            node2 = [node2; pos(iy(nc),:)];
            strng = [strng; iv(nc)];
        end
    end
end

end

function data = drawnodes(data, N)
% Node plotter. N = (90,1) with 1s for nodes to plot and 0s to ignore.
%
% 
hold on;
pos = data.sourcemodel.pos;
%v   = pos*0.9;

if isfield(data.sourcemodel,'vi')
    % if pos is cell, it's because we've passed both a vertex set and a
    % vector that describes which vertex belongs to which roi
    v  = pos;
    vi = data.sourcemodel.vi;
    n  = unique(vi);
    for i = 1:length(n)
        these = find(vi==n(i));
        roi(i,:) = spherefit(v(these,:));
    end
    pos = roi;
    for ip = 1:length(pos)
        [~,this]  = min(cdist(pos(ip,:),data.mesh.vertices));
        pos(ip,:) = data.mesh.vertices(this,:);
    end
end


bounds = [min(data.mesh.vertices); max(data.mesh.vertices)];
offset = 0.99;
for ip = 1:3
    pos(:,ip) = bounds(1,ip) + ((bounds(2,ip)-bounds(1,ip))) .* ...
                (pos(:,ip) - min(pos(:,ip)))./(max(pos(:,ip)) - min(pos(:,ip)));
    pos(:,ip) = pos(:,ip)*offset;
end

% redirect to clseast mesh point (vertex?)
for ip = 1:length(pos)
    [~,this]  = min(cdist(pos(ip,:),data.mesh.vertices));
    pos(ip,:) = data.mesh.vertices(this,:);
end

v = pos;

if size(N,1) > 1 && size(N,2) > 1
    cols = {'r' 'm','y','g','c','b'};
    if size(size(N,2)) == 90
        N = N';
    end
    
    for j = 1:size(N,2)
        ForPlot = v(find(N(:,j)),:) + (1e-2 * (2*j) ) ;
        s       = find(N);
        col     = cols{j};
        for i   = 1:length(ForPlot)
            scatter3(ForPlot(i,1),ForPlot(i,2),ForPlot(i,3),70,col,'filled',...
                'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6);        hold on;
        end
    end
    
elseif iscell(N)
    % pass a cell of length sm with color strings
    % e.g. N=cell(90,1);
    %      N{1} = 'r';
    ForPlot = [];
    for j = 1:length(N)
        if any(N{j})
            v0 = v(j,:);
            ForPlot = [ForPlot v0];
            scatter3(v0(1),v0(2),v0(3),150,N{j},'filled');
        end
    end
    
else
    ForPlot = v(find(N),:);
    %s       = find(N);
    s = ones(length(find(N)),1)*150;
    for i   = 1:size(ForPlot,1)
        col = 'r';
        scatter3(ForPlot(i,1),ForPlot(i,2),ForPlot(i,3),s(i),'r','filled');
    end
end
%RGB = makecolbar(ForPlot);
%set(gcf,'DefaultAxesColorOrder',RGB); jet;
colorbar

data.drawnodes.data = ForPlot;

end

function RGB = makecolbar(I,netcmap)
% Register colorbar values to an overlay /  T-vector
%

if any(any(netcmap ~= 0))
    Colors = colormap(netcmap);
else
    Colors   = jet;
end

NoColors = length(Colors);

Ireduced = (I-min(I))/(max(I)-min(I))*(NoColors-1)+1;
RGB      = interp1(1:NoColors,Colors,Ireduced);

end


function [y,data] = parse_overlay(x,data)
% Functinal overlay data parser: volumes, filenames, statements
%
%
%

if iscell(x)
    % multiple overlays (e.g. {'curvature','nif.nii'} )
   for i = 1:length(x)
       [y{i},data] = parse_overlay(x{i},data);
   end
   return;
end

if ischar(x)
    [fp,fn,fe] = fileparts(x);
    
    
    switch fn
        case 'curvature'
            % since we've precomputed curvature anyway, allow calling 
            % ['overlay', 'curvature'] & return it as y
            y = data.mesh.curvature;
            return;
    end
       
    switch fe
        
        case{'.gz'}
            fprintf('Unpacking .gz\n');
            gunzip(x);
            x = strrep(x,'.gz','');
            
            [y,data] = parse_overlay(x,data);
            return;
            
        case{'.nii'}
            
            % add waitbar
            wb = waitbar(0,'Reading volume: Please wait...');
            
            % load nifti volume file
            fprintf('Reading Nifti volume\n');
%             ni    = load_nii(x);
%             vol   = ni.img;

            ni  = ft_read_mri(x);
            vol = ni.anatomy; 

            R = ni.transform(1:3,1:3);
            t = ni.transform(1:3,end);

            % real-world cooredinates = 
            %[x; y; z] = R * [i; j; k] + t

            [nx,ny,nz] = size(vol);

            % compute edges
            E1 = R*[1 ;1 ;1 ]+t;
            E2 = R*[nx;ny;nz]+t;

            % real-world cartesian coordinates
            xarr = linspace(E1(1),E2(1),nx);
            yarr = linspace(E1(2),E2(2),ny);
            zarr = linspace(E1(3),E2(3),nz);
            
            % retain header info?
            data.volume.fname     = x;
            data.volume.hdr       = ni.hdr;
            data.volume.vol       = vol;
            
            % NEW: subsample volume by 2. This speeds up subsequent patch
            % reduction enormously without at relatively little cost to
            % resolution
            [xg,yg,zg] = meshgrid(yarr,xarr,zarr);
            [NX, NY, NZ, vol] = reducevolume(xg,yg,zg,vol,2);
            xarr = linspace(xarr(1),xarr(end),size(vol,1));
            yarr = linspace(yarr(1),yarr(end),size(vol,2));
            zarr = linspace(zarr(1),zarr(end),size(vol,3));

            % retain coordinate system to pass to vol2surf
            data.volume.grid.x = xarr;
            data.volume.grid.y = yarr;
            data.volume.grid.z = zarr;
            
            
            [y,data] = vol2surf(vol,data,wb);
            
            
            % ensure sourcemodel (pos) is around same scale as mesh boundaries
            m = min(data.mesh.vertices);% *1.1;
            M = max(data.mesh.vertices);% *1.1;

            pos      = data.sourcemodel.pos;
            V        = pos - repmat(spherefit(pos),[size(pos,1),1]);
            V(:,1)   = m(1) + ((M(1)-m(1))).*(V(:,1) - min(V(:,1)))./(max(V(:,1)) - min(V(:,1)));
            V(:,2)   = m(2) + ((M(2)-m(2))).*(V(:,2) - min(V(:,2)))./(max(V(:,2)) - min(V(:,2)));
            V(:,3)   = m(3) + ((M(3)-m(3))).*(V(:,3) - min(V(:,3)))./(max(V(:,3)) - min(V(:,3)));
            data.sourcemodel.pos = V;
            
        case{'.gii'}
            % load gifti functional
            gi = gifti(x);
            y  = double(gi.cdata);
            if length(y) ~= length(data.sourcemodel.pos)
                fprintf('Gifti overlay does not match sourcemodel!\n');
                if length(y) == length(data.mesh.vertices)
                    fprintf(['...but it does match the mesh size.\nUsing '...
                    'mesh vertices as sourcemodel\n']);
                    data.sourcemodel.pos = data.mesh.vertices;
                end     
            end
    end
end

if isnumeric(x) && ndims(x)==3
    % this is a pre-loaded nifti volume
    fprintf('This is a pre-loaded 3D nifti volume: extracting...\n');
    
    % add waitbar
    wb = waitbar(0,'Preloaded volume: Please wait...');
    
    % NEW: subsample volume by 2 to speed up patch reduction
    fprintf('Subsampling structural volume\n');
    [NX, NY, NZ, x] = reducevolume(x,2);

    [y,data] = vol2surf(x,data,wb);
    
end

end

function y = sym_pad_vector(x,n)
% Symmetrically pad a vector with zeros

if length(x) ~= n
    k = n - length(x);
    k = floor(k/2);
    y = [zeros(1,k) x zeros(1,k)];
    
else y = x;
end

end

function [y,data] = vol2surf(vol,data,wb)
% Convert a FUNCTIONAL volume to surface. Delauay triangulate, compute a
% set of 'source' vertices. Take liberties with downsampling and smoothing


% bounds:
S = size(vol);

% check if it's a 'full' volume!
if length(find(vol)) == prod(S)
    vol = vol - mode(vol(:));
end

% waitbar
waitbar(.15,wb,'Reading volume: Smoothing volume');

% a little smoothing
vol = smooth3(vol,'gaussian');

% New --- only exists is vol was loaded by atemplate
try
    pixdim = data.volume.hdr.dime.pixdim(2:4);
end

if ~isfield(data,'volume') || ~isfield(data.volume,'grid')
    fprintf('Using computed voxel coordinates\n');
    x = 1:size(vol,1);
    y = 1:size(vol,2);
    z = 1:size(vol,3);
elseif isfield(data.volume,'grid')
    fprintf('Using real-world voxel coordintes from nifti\n');
    x = data.volume.grid.x;
    y = data.volume.grid.y;
    z = data.volume.grid.z;
end

% waitbar
waitbar(.30,wb,'Reading volume: Indexing volume');

% find indices of tissue in old grid
[nix,niy,niz] = ind2sub(size(vol),find(vol));
[~,~,C]       = find(vol);

% waitbar
waitbar(.45,wb,'Reading volume: Compiling new vertex list');

% compile a new vertex list
fprintf('Compiling new vertex list (%d verts)\n',length(nix));
v = [x(nix); y(niy); z(niz)]';
v = double(v);
try
    v = v*diag(pixdim);
end

% apply affine if req.
if isfield(data.overlay,'affine')
    affine = data.overlay.affine;
    if length(affine) == 4
        fprintf('Applying affine transform\n');
        va = [v ones(length(v),1)]*affine;
        v  = va(:,1:3);
    end
end

% Fit this gridded-volume inside the box extremes of the mesh
B        = [min(data.mesh.vertices); max(data.mesh.vertices)];
V        = v - repmat(spherefit(v),[size(v,1),1]);
V(:,1)   = B(1,1) + ((B(2,1)-B(1,1))).*(V(:,1) - min(V(:,1)))./(max(V(:,1)) - min(V(:,1)));
V(:,2)   = B(1,2) + ((B(2,2)-B(1,2))).*(V(:,2) - min(V(:,2)))./(max(V(:,2)) - min(V(:,2)));
V(:,3)   = B(1,3) + ((B(2,3)-B(1,3))).*(V(:,3) - min(V(:,3)))./(max(V(:,3)) - min(V(:,3)));
v        = V;


% % new grid
% fprintf('Generating grid for volume data\n');
% x = linspace(B(1,1),B(2,1),S(1));
% y = linspace(B(1,2),B(2,2),S(2));
% z = linspace(B(1,3),B(2,3),S(3));
% 
% % find indiced of tissue in old grid
% [nix,niy,niz] = ind2sub(size(vol),find(vol));
% [~,~,C]       = find(vol);
% 
% % compile a new vertex list
% fprintf('Compiling new vertex list (%d verts)\n',length(nix));
% v = [x(nix); y(niy); z(niz)]';
% v = double(v);


% reduce patch
fprintf('Reducing patch density\n');

% waitbar
waitbar(.50,wb,'Reading volume: Triangulating');

nv  = length(v);
tri = delaunay(v(:,1),v(:,2),v(:,3));
fv  = struct('faces',tri,'vertices',v);
count  = 0;

% waitbar
waitbar(.60,wb,'Reading volume: Smoothing');

% smooth overlay at triangulated points first
Cbound = [min(C) max(C)];
C      = spm_mesh_smooth(fv,double(C),4);
C      = Cbound(1) + (Cbound(2)-Cbound(1)).*(C - min(C))./(max(C) - min(C));

% waitbar
waitbar(.70,wb,'Reading volume: Reducing patch density');

while nv > 10000
   fv  = reducepatch(fv, 0.5);
   nv  = length(fv.vertices);
   count = count + 1;
end

% print
fprintf('Patch reduction finished\n');
fprintf('Using nifti volume as sourcemodel and overlay!\n');
fprintf('New sourcemodel has %d vertices\n',nv);

% waitbar
waitbar(.90,wb,'Reading volume: Computing colours for reduced patch');

% find the indices of the retained vertexes only
fprintf('Retrieving vertex colours\n');
Ci = compute_reduced_indices(v, fv.vertices);

% waitbar
waitbar(1,wb,'Reading volume: eComplete');
close(wb);


% Update sourcemodel and ovelray data
v                    = fv.vertices;
data.sourcemodel.pos = v;
y                    = C(Ci);

end


function indices = compute_reduced_indices(before, after)
% Compute the retained functional (colour) values after patch reduction
%

indices = zeros(length(after), 1);
for i = 1:length(after)
    dotprods = (before * after(i, :)') ./ sqrt(sum(before.^2, 2));
    [~, indices(i)] = max(dotprods);
end
end

function y = rescale(x,S)

y = S(1) + (S(2)-S(1)) .* (x - min(x) ) / ...
    ( max(x) - min(x) );

end

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
   [L,data] = parse_overlay(L,data);
   
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
           [data,i]   = sort_template(data,i);
           L          = i.L;
       end
   end
end

% method for searching between the 3D coordinate systems
%-------------------------------------------------------------
if ischar(data.overlay.method)
    if ismember(data.overlay.method,{'euclidean','spheres','precomputed (AAL)','raycast','aal','aal_light'})
         method = data.overlay.method;
    else,method = 'euclidean';  
    end
else
    method = data.overlay.method{1};
end

% If requested multi-overlay, where one is curvature, draw the curvature
%--------------------------------------------------------------------------
if iscell(L)
    % curvature is first
    C = L{1};
        
    % spm mesh smoothing
    fprintf('Smoothing overlay...\n');
    y = spm_mesh_smooth(data.mesh, double(C(:)), 4);
    percNaN = length(find(isnan(C)))/length(C)*100;
    newpNaN = length(find(isnan(y)))/length(y)*100;
    
    % when using a NaN-masked overlay, smoothing can result in all(nan) or
    % an increase in the number of nans: enforce a 5% tolerance on this, which
    % forces reverting to the uns-smoothed version if reached
    if all(isnan(y)) || newpNaN > (percNaN*1.05)
        fprintf('Reverting to non-smoothed overlay due to too many NaNs\n');
        y = C(:);
    end
    
    set(data.mesh.h,'FaceVertexCData',y(:),'FaceColor','interp');
    
    drawnow;
    shading interp
    % force symmetric caxis bounds
    s = max(abs(y(:))); caxis([-s s]);
    fprintf('Generating combined function+curvature colormap\n');
    themap = [flipud(gray(256)); flipud(jet(256))];
    alpha 1;

    % replace L and flag requirement for a new axis
    L = L{2};
    NewAx = 1;
elseif isfield(data.overlay,'NewAx')
    NewAx = data.overlay.NewAx;
else
    NewAx = 0;
end

data.overlay.NewAx = NewAx; % save this in case of recursive calls

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

if NewAx
    % if overlaying functional and curvature
    % *no longer generates new axes!*
    ax1     =  gca;
    ax1_pos = get(ax1,'Position'); 
    %ax2 = axes('Position',ax1_pos,'visible',0);  
    h0  = findobj(ax1,'type','patch');
end
    
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

% first-pass attempt at ensuring alignment of sourcemodel and outer mesh
%pos = fixmesh(mesh,pos);

% Overlay
v  = pos;                       % sourcemodel vertices
x  = v(:,1);                    % AAL x verts
mv = mesh.vertices;             % brain mesh vertices
nv = length(mv);                % number of brain vertices
S  = [min(L(:)),max(L(:))];     % min max values
S0 = max(abs(L(:)));

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
        UseFaceVertex = 'vertex'; 
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
                vol(dv(i,1)-ndv(1),dv(i,2)-ndv(2),dv(i,3)-ndv(3)) = ...
                vol(dv(i,1)-ndv(1),dv(i,2)-ndv(2),dv(i,3)-ndv(3)) + a(I);                
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
                %fcol = [fcol; S(1); S(2)];
                fcol  = S(1) + ((S(2)-S(1))).*(fcol - min(fcol))./(max(fcol) - min(fcol));
                %fcol = fcol(1:end-2);
                
                if ~NewAx
                    set(mesh.h,'FaceVertexCData',fcol(:),'FaceColor','interp');
                else
                    % if we're overlaying the functional ray-casted colours on top of the curvature
                    if ~isempty(data.overlay.thresh)
                          thrsh = data.overlay.thresh;
                    else; thrsh = .4;
                    end
                    
                    fcol  = S(1) + ((S(2)-S(1))).*(fcol - min(fcol))./(max(fcol) - min(fcol));
                    
                    fcol_orig = fcol;
                    thr  = max(abs(fcol))*thrsh;
                    inan = find(abs(fcol) < thr);
                    falpha = 1*ones(size(fcol));
                    falpha(inan)=0;
                    %fcol = fcol.*falpha;
                                        
                    fprintf('Rescaling overlay values\n');
                    
                    % get curvature colours
                    y0 = data.mesh.h.FaceVertexCData;
                    
                    % rescale y0 into colour map part 1:
                    % 1:256
                    y0 = 256 * (y0-min(y))./(max(y0)-min(y0));
                    
                    % rescale fcol into colour map part 2:
                    % 257:512
                    m = 257;
                    n = 512;
                    fcol = [fcol; -S0; S0];
                    fcol = n + (m - n) .* (fcol-min(fcol))./(max(fcol)-min(fcol));
                    fcol = fcol(1:end-2);
                    
                    % mask new functional colours
                    fcol = fcol.*falpha;
                    
                    new_over = y0;
                    these    = find(fcol);
                    new_over(these) = fcol(these); 
                    
                    set(data.mesh.h,'FaceVertexCData',new_over(:),'FaceColor','interp');
                    colormap(themap);
                    data.overlay.themap=themap;
                end
        end
        
        % Use symmetric colourbar and jet as defaults
        if NewAx
            %colormap(ax2,'jet');
            %s = max(abs(fcol(:))); caxis(ax2,[-s s]);
            caxis([0 512]);
        else
            colormap('jet');
            s = max(abs(fcol(:))); caxis([-s s]);
            alpha 1;
        end
        %alpha 1;
        
        % Return the face colours
        data.overlay.data  = fcol(:);       % the functional vector
        data.overlay.steps = step;          % the depths at which searched
        data.overlay.hits  = hits;          % num hits / intersects at each depth
        data.overlay.cast  = UseFaceVertex; % whether computed for faces or vertices
        
        if NewAx
            data.overlay.data_nonan = fcol_orig;
        end
        
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

        % for no interpolation, set w = 1;
        % w = 1;
        
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
        
        
    case {'aal','aal90','aal_90'}
        % project into pre-computed AAL parcellation - 1 value per region
        %
        
        % new method for AAL90: wrapper on ray casting routine
        load DenseAAL.mat

        ol    = zeros(length(v),1);
        for i = 1:length(L)
            these = find(vi==i);
            ol(these) = L(i);
        end
        
        % update sourcemodel
        v = fit_check_source2mesh(v,data.mesh); 
        data.sourcemodel.pos = v;
        data.overlay.orig    = ol;
        data.overlay.method  = data.overlay.method{2};
        
        data = overlay(data,ol,write,fname,colbar);
        return;
        
    case {'aal_light','aal_reduced','aal_red'};
        % project into pre-computed AAL parcellation - 1 value per region
        % - reduce dversion
        
        % new method for AAL90: wrapper on ray casting routine
        load LightAAL.mat

        ol    = zeros(length(v),1);
        for i = 1:length(L)
            these = find(vi==i);
            ol(these) = L(i);
        end
        
        % update sourcemodel
        v = fit_check_source2mesh(v,data.mesh); 
        data.sourcemodel.pos = v;
        data.overlay.orig    = ol;
        data.overlay.method  = data.overlay.method{2};
        data.overlay.atlasvalues = L;
        data = overlay(data,ol,write,fname,colbar);
        return;
        
        
    case {'harvox','ho','harvardoxford','hoa','harvard_oxford'}
        % project into pre-computed Harvard-Oxford Atlas
        % parcellation - 1 value per region
        % then use interp/project method above to render on user mesh
        
        % new method for AAL90: wrapper on above routine
        load HarvOx.mat

        ol    = zeros(length(v),1);
        for i = 1:length(L)
            these = find(vi==i);
            ol(these) = L(i);
        end
        
        % update sourcemodel
        v = fit_check_source2mesh(v,data.mesh); 
        data.sourcemodel.pos = v;
        data.overlay.orig    = ol;
        data.overlay.method  = data.overlay.method{2};
        
%         % compute roi centres for labelling if req
%         data.overlay.atlas_flag = 1;
%         %v_roi = get_roi_centres0(v,vi);
%         data.overlay.atlas.labels = labels(:,2);
%         %data.overlay.rois = v_roi;
%         data.overlay.atlas.v  = v;
%         data.overlay.atlas.vi = vi;
        
        %v_roi(isnan(v_roi)) = 0;
        
        data = overlay(data,ol,write,fname,colbar);
        return;
        
    case {'aal116','aal_116'};
        % project into pre-computed AAL parcellation - 1 value per region
        %
        
        load AAL116.mat

        ol    = zeros(length(v),1);
        for i = 1:length(L)
            these = find(vi==i);
            ol(these) = L(i);
        end
        
        % update sourcemodel
        v = fit_check_source2mesh(v,data.mesh); 
        data.sourcemodel.pos = v;
        data.overlay.orig    = ol;
        data.overlay.method  = data.overlay.method{2};
        
        data = overlay(data,ol,write,fname,colbar);
        return;
        
    case {'hoa_cerebellum','hoac','hoa_c'};
        % project into pre-computed AAL parcellation - 1 value per region
        %
        
        load HOA_Cerebellum.mat

        ol    = zeros(length(v),1);
        for i = 1:length(L)
            these = find(vi==i);
            ol(these) = L(i);
        end
        
        % update sourcemodel
        v = fit_check_source2mesh(v,data.mesh); 
        data.sourcemodel.pos = v;
        data.overlay.orig    = ol;
        data.overlay.method  = data.overlay.method{2};
        
        data = overlay(data,ol,write,fname,colbar);
        return;
        
        
    case {'user'}
        % project into user specified parcellation - 1 value per region
        %
        
        v  = data.sourcemodel.pos;
        vi = data.sourcemodel.vi; 

        ol    = zeros(length(v),1);
        for i = 1:length(L)
            these = find(vi==i);
            ol(these) = L(i);
        end
        
        % update sourcemodel
        %v = fit_check_source2mesh(v,data.mesh); 
        data.sourcemodel.pos = v;
        data.overlay.orig    = ol;
        data.overlay.method  = data.overlay.method{2};
        
        data = overlay(data,ol,write,fname,colbar);
        return;
end


switch method
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
        %------------------------------------------------------------------
        fprintf('Smoothing overlay...\n');
        y  = spm_mesh_smooth(mesh, y(:), 4);
        y(isnan(y)) = 0;
        y  = S(1) + ((S(2)-S(1))).*(OL - min(OL))./(max(OL) - min(OL));
        y(isnan(y)) = 0;

        % return these in data structre
        try data.overlay.data           = y; end
        try data.overlay.smooth_weights = M; end
        try data.overlay.NumComp        = NumComp; end
        try data.overlay.indz           = indz;    end
        try data.overlay.w              = w;       end

        if ~NewAx
            set(mesh.h,'FaceVertexCData',y(:),'FaceColor','interp');
        else
            
            if ~isempty(data.overlay.thresh)
                thrsh = data.overlay.thresh;
            else; thrsh = .4;
            end
            fcol      = y;
            fcol  = S(1) + ((S(2)-S(1))).*(fcol - min(fcol))./(max(fcol) - min(fcol));
                    
            fcol_orig = fcol;
            thr  = max(abs(fcol))*thrsh;
            inan = find(abs(fcol) < thr);
            falpha = 1*ones(size(fcol));
            falpha(inan)=0;
            %fcol = fcol.*falpha;
            
            fprintf('Rescaling overlay values\n');
            
            % get curvature colours
            y0 = data.mesh.h.FaceVertexCData;
            
            % rescale y0 into colour map part 1:
            % 1:256
            y0 = 256 * (y0-min(y))./(max(y0)-min(y0));
            
            % rescale fcol into colour map part 2:
            % 257:512
            m = 257;
            n = 512;
            fcol = [fcol; -S0; S0];
            fcol = n + (m - n) .* (fcol-min(fcol))./(max(fcol)-min(fcol));
            fcol = fcol(1:end-2);
            
            % mask new functional colours
            fcol = fcol.*falpha;
            
            new_over = y0;
            these    = find(fcol);
            new_over(these) = fcol(these);
            
            set(data.mesh.h,'FaceVertexCData',new_over(:),'FaceColor','interp');
            colormap(themap);
            data.overlay.themap=themap;
        end
        drawnow;
        shading interp
        % force symmetric caxis bounds
        s = max(abs(y(:))); caxis([-s s]);
        if ~NewAx
            colormap('jet');
        else
            %colormap(ax2,'jet');
            caxis([0 512]);
        end
        alpha 1;
        
end

S = [-s s];

% post-hoc template reduction here?
%--------------------------------------------------------------------------
if isfield(data,'post_parcel')
    fprintf('Also computing requested (parcellated) atlas transform\n');
    newv = data.post_parcel{1};
    newv = fit_check_source2mesh(newv,struct('vertices',data.sourcemodel.pos));
    D0   = cdist(data.sourcemodel.pos,newv);
    
    % assume the vertex value of the closest mesh point
    for i = 1:size(D0,2)
        [~,ind] = min(D0(:,i));
        newL(i) = data.overlay.orig(ind);
    end
    newL = [newL S(1) S(2)];
    newL = S(1) + ((S(2)-S(1))).*(newL - min(newL))./(max(newL) - min(newL));
    newL = newL(1:end-2);
    
    % if the second input was supplied, compute average atlas value as mean
    % of parcel vertices
    try 
        iv = data.post_parcel{2};
        fprintf('Computing parcel means\n');
        n = unique(iv);
        n(n==0) = [];
        ParVal  = iv*0; 
        for i = 1:length(n)
            these = find(iv==n(i));
            ParcelMean(i) = mean( newL(these) );
            MeanLoc(i,:)  = mean([spherefit( newv(these,:) ) ; mean( newv(these,:) )]);
            ParVal(these) = ParcelMean(i);
        end
        ParcelMean = S(1) + ((S(2)-S(1))).*(ParcelMean - min(ParcelMean))./(max(ParcelMean) - min(ParcelMean));
        ParVal     = S(1) + ((S(2)-S(1))).*(ParVal - min(ParVal))./(max(ParVal) - min(ParVal));
    end
    
    % box bound new set
    %MeanLoc = fit_check_source2mesh(MeanLoc,data.mesh);
    
    % return the new v and data
    data.post_parcel      = [];
    data.post_parcel.pos  = newv;
    data.post_parcel.data = newL;  

    try 
        % return atlas data as means
        data.post_parcel.ParcelMean = ParcelMean;
        data.post_parcel.ParcelCent = MeanLoc;
        data.post_parcel.ParVal     = ParVal;
        data.post_parcel.ParcelID   = iv;
    end 
end



if colbar && ~NewAx
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
    V = genvol(data.mesh.vertices,data.overlay.data,[256 256 256]);
    %new = mesh;
    %dim = ceil(nthroot(length(y),3));
    %V   = sm2vol(new.vertices,dim*3,y,256);
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

function pos = fit_check_source2mesh(pos,g)
% check box bounds of sourcemodel and mesh are matched
%
% pos = nx3 source model, g = mesh gifti/struct

% ensure sourcemodel (pos) is around same scale as mesh boundaries
m = min(g.vertices);% *1.1;
M = max(g.vertices);% *1.1;

% refit box boundaries
V        = pos - repmat(spherefit(pos),[size(pos,1),1]);
V(:,1)   = m(1) + ((M(1)-m(1))).*(V(:,1) - min(V(:,1)))./(max(V(:,1)) - min(V(:,1)));
V(:,2)   = m(2) + ((M(2)-m(2))).*(V(:,2) - min(V(:,2)))./(max(V(:,2)) - min(V(:,2)));
V(:,3)   = m(3) + ((M(3)-m(3))).*(V(:,3) - min(V(:,3)))./(max(V(:,3)) - min(V(:,3)));
pos      = V;

end

function y = makeodd(x)
% nearest odd number
%
%

y = 2*floor(x/2)+1;

end


function x = killinterhems(x);
% kill interhemispherics, assuming a square adjacency matrix
%
%

S  = size(x);
xb = (S(1)/2)+1:S(1);
yb = (S(2)/2)+1:S(2);
xa = 1:S(1)/2;
ya = 1:S(2)/2;

x(xa,yb) = 0;
x(xb,ya) = 0;

end




function newpos = fixmesh(g,pos)
% redirect sourcemodel points onto closest mesh vertex
%
% AS

v = g.vertices;
v = v - repmat(spherefit(v),[size(v,1),1]); % Centre on ~0
g.vertices=v;

% Centre on ~0
pos = pos - repmat(spherefit(pos),[size(pos,1),1]);

newpos = pos;
for i = 1:length(pos)
    this  = pos(i,:);
    [t,I] = maxpoints(cdist(v,this),1,'max');
    newpos(i,:) = v(I,:);    
end

end

function [g,pos,h,p] = meshmesh(g,write,fname,fighnd,a,pos,hemisphere,affine,flip,inflate,checkori,dofillholes,optimise)
% Main brain-mesh plot. 
%
% Separates hemispheres, computes curvature, check orientation, inflate, 
% centre & fill holes. Plots the patch as a semi-transparent glass brain 
% for other plot functiopns to add to. 
%
% Returns everything in data.mesh
%

if iscell(pos)
   orig_pos = pos;
   pos      = pos{1};
end

if isempty(a);
    a = .6;
end

v = g.vertices;

% apply affine if req.
if length(affine) == 4
    fprintf('Applying affine transform\n');
    va = [v ones(length(v),1)]*affine;
    v  = va(:,1:3);
    g.vertices = v;
end

% flip x/y if required but POST affine transform
if flip
    v = g.vertices;
    v = v(:,[2 1 3]);
    g.vertices = v;
end


% check rotation
yl = max(v(:,2)) - min(v(:,2));
xl = max(v(:,1)) - min(v(:,1));

if xl > yl
    v       = v(:,[2 1 3]);
    g.faces = g.faces(:,[2 1 3]);
end

% centre and scale mesh
g.vertices = v - repmat(spherefit(v),[size(v,1),1]);


% ensure sourcemodel (pos) is around same scale as mesh boundaries
pos = fit_check_source2mesh(pos,g);

% explicitly optimise the alignment of the cloud points?
% (source model vertices and mesh vertices)
if optimise
    pos = align_clouds_3d(g.vertices,pos);
    clf; drawnow;
end

% m = min(g.vertices);% *1.1;
% M = max(g.vertices);% *1.1;
% 
% V        = pos - repmat(spherefit(pos),[size(pos,1),1]);
% V(:,1)   = m(1) + ((M(1)-m(1))).*(V(:,1) - min(V(:,1)))./(max(V(:,1)) - min(V(:,1)));
% V(:,2)   = m(2) + ((M(2)-m(2))).*(V(:,2) - min(V(:,2)))./(max(V(:,2)) - min(V(:,2)));
% V(:,3)   = m(3) + ((M(3)-m(3))).*(V(:,3) - min(V(:,3)))./(max(V(:,3)) - min(V(:,3)));
% pos      = V;

% % calculate curvature for shading
curv = docurvature(struct('vertices',g.vertices,'faces',g.faces));

% inflate
if inflate
    fprintf('Inflating mesh\n');
    g = spm_mesh_inflate(struct('vertices',g.vertices,'faces',g.faces),150);
end


if checkori
    b = CheckOrientationMesh(g);
    waitfor(b)
    g = b;
end


% only one hemisphere?
v = g.vertices;
f = g.faces;
c = spherefit(v);

% save centre point for later - e.g. label setting
g.centre = c;

left  = find(v(:,1) < c(1));
right = find(v(:,1) > c(1));

lfaces = find(sum(ismember(f,left),2)==3);
rfaces = find(sum(ismember(f,right),2)==3);

% return left/right indices
g.vleft            = v*NaN;
g.vleft(left,:)    = v(left,:);
g.vright           = v*NaN;
g.vright(right,:)  = v(right,:);
g.fleft            = f*NaN;
g.fleft(lfaces,:)  = f(lfaces,:);
g.fright           = f*NaN;
g.fright(rfaces,:) = f(rfaces,:);
g.curvature        = curv;



% fill holes stemming from hemisphere separation - if requested
%---------------------------------------------------------------
if dofillholes
    fprintf('Filling holes resulting from hemisphere separation..\n');
    
    % left side
    x = [];
    holecellarray = findTriMeshHoles(g.fleft,g.vleft);
    for i = 1:length(holecellarray)
        edges   = holecellarray{i};  % indices of vertices around edge
        ek      = nchoosek(edges,3); % triangles at all perimiter points
        x       = [x; ek];           % new faces
    end
    
    x0 = x;
    
    % update faces
    g.fleft = [g.fleft; x];
        
    % right side
    x = [];
    holecellarray = findTriMeshHoles(g.fright,g.vright);
    for i = 1:length(holecellarray)
        edges   = holecellarray{i};  % indices of vertices around edge
        ek      = nchoosek(edges,3); % triangles at all perimiter points
        x       = [x; ek];           % new faces
    end
    
    % update faces
    g.fright = [g.fright; x];
    
    % update full-model faces
    g.faces = [g.faces; x0; x];
end


switch hemisphere
    case{'left','L','l'}
        pg.vertices = g.vleft;
        pg.faces    = g.fleft;
        
        %pg.vertices         = v*NaN;
        %pg.vertices(left,:) = v(left,:);
        %pg.faces           = f*NaN;
        %pg.faces(lfaces,:) = f(lfaces,:);
  
    case{'right','R','r'}
        pg.vertices = g.vright;
        pg.faces    = g.fright;
        
        %pg.vertices          = v*NaN;
        %pg.vertices(right,:) = v(right,:);
        %pg.faces           = f*NaN;
        %pg.faces(rfaces,:) = f(rfaces,:);        
        
    otherwise
        pg = g;
end

% sanity check
if any( min(pg.faces(:)) == 0)
    bad = find(pg.faces == 0);
    pg.faces(bad) = NaN;
end

pg.faces = double(pg.faces);
pg.faces(pg.faces==0) = NaN;

% plot
if ~isempty(fighnd)
    if isnumeric(fighnd)
        % Old-type numeric axes handle
        h = plot(fighnd,gifti(pg));
    elseif ishandle(fighnd)
        % new for matlab2017b etc
        % [note editted gifti plot function]
        h = plot(gifti(pg),'fighnd',fighnd);
    end
else
    h  = plot(gifti(pg));
end
C = [.5 .5 .5];

set(h,'FaceColor',[C]); box off;
grid off;  set(h,'EdgeColor','none');
alpha(a); set(gca,'visible','off');

%h = get(gcf,'Children');
%set(h(end),'visible','off');
set(gca,'visible','off')

drawnow; hold on;

p = [];

if write == 1
    fprintf('Writing mesh gifti file: %s\n',[fname '.gii']);
    gout.vertices = g.vertices;
    gout.faces    = g.faces;
    gout = gifti(gout);
    save(gout,fname);
end



if exist('orig_pos')
    % is pos was a cell, re-instate 
    orig_pos{1} = pos;
    pos         = orig_pos;
end

end

function C = docurvature(M)
% Computes the curvature of a mesh
%
%

A = spm_mesh_adjacency(M);
A = sparse(1:size(M.vertices,1),1:size(M.vertices,1),1./sum(A,2)) * A;

C = (A-speye(size(A))) * double(M.vertices);
N = spm_mesh_normals(M);
C = sign(sum(N.*C,2)) .* sqrt(sum(C.*C,2));

end

function data = addlabels(data,V,all_roi_tissueindex,thelabels)
% Add labels to the plot.
%
% If using AAL90 sourcemodle, these are automatic.
%
% If using another sourcemodel:
% - provide the all_roi_tissueindex from fieldtirp. This is a
% 1xnum_vertices vector containing indices of rois (i,e. which verts belong
% to which rois).
% Also provide labels!
%
pos = data.sourcemodel.pos;

if ( ~isempty(thelabels) && ~isempty(all_roi_tissueindex) ) &&...
   ( length(pos) == length(all_roi_tissueindex) ) &&...
   ( length(thelabels) == length(unique(all_roi_tissueindex(all_roi_tissueindex~=0))) )
    
    labels = strrep(thelabels,'_',' ');
    v      = get_roi_centres(pos,all_roi_tissueindex);
    roi    = all_roi_tissueindex;
    
elseif length(V) == 90
    
    load('AAL_labels');
    labels = strrep(labels,'_',' ');
    v      = pos*0.95;
    roi    = 1:90;
elseif (length(V) == length(thelabels)) &&...
       (length(V) == length(pos))
    
   labels = strrep(thelabels,'_',' ');
    v = pos*0.95;
    roi = 1:length(V);

% elseif isfield(data,'overlay') && isfield(data.overlay,'atlas_flag') 
%     % user is calling one of the built-in atlas parcellations, for which we
%     % already know the labels and label locations... (AAL or HOA)
%     
%     v  = data.overlay.atlas.v;
%     vi = data.overlay.atlas.vi;
%     labels = data.overlay.atlas.labels;
%     
%     fprintf('Translating label positions to centres of parcellated ROIs\n');
%     fprintf('(this can take a minute)\n');
%     label_rois = diag(vi)*data.overlay.smooth_weights;
%     
%     for i = 1:length(vi)
%         roin = find(vi==i);
%         [x,y]=find(label_rois(roin,:));
%         these(unique(y)) = i;
%     end
%     
%     % now compute parcel centres for each roi
%     rois = get_roi_centres(data.mesh.vertices,these);
%     v = rois;
%     roi = 1:length(v);
%     n = length(v);
%     V = ones(n,n);
%     v = double(v);
        
else
    fprintf('Labels info not right!\n');
    return
end

data.labels.roi     = roi;
data.labels.labels  = labels;
data.labels.centres = v;

% compile list of in-use node indices
%--------------------------------------------------------------------------
to = []; from = []; V(isnan(V)) = 0;
for i  = 1:size(V,1)
    ni = find(logical(V(i,:)));
    if any(ni)
        to   = [to   roi(ni)];
        from = [from roi(repmat(i,[1,length(ni)])) ];
    end
end

AN  = unique([to,from]);
AN  = AN(AN~=0);
off = 1.5;
data.labels.in_use = AN;

if sum(v(:,3)) == 0
     off3 = 0;
else;off3 = off;
end

% recompute brain centre
c = spherefit(data.mesh.vertices);

% add these to plot with offset
%--------------------------------------------------------------------------
for i = 1:length(AN)
    L = labels{AN(i)};
    pos = v(AN(i),:);
    
    if pos(1) < c(1)
        hemi = 'L';
    elseif pos(1) > c(1)
        hemi = 'R';
    else
        hemi = ' ';
    end
    
    switch hemi
        case 'L';
            t(i) = text(v(AN(i),1)-(off*5),v(AN(i),2)-(off*5),v(AN(i),3)+off3,L);
        case 'R';
            t(i) = text(v(AN(i),1)+(off*2),+v(AN(i),2)+(off*2),v(AN(i),3)+off3,L);
        otherwise
            t(i) = text(v(AN(i),1),v(AN(i),2),v(AN(i),3),L);
    end
end
try
    % this fails if the network was empty!
    set(t,'Fontsize',14)
end

end

function [C,verts] = get_roi_centres(pos,all_roi_tissueindex)
% Find centre points of rois
%
%
v   = pos;
roi = all_roi_tissueindex;

i   = unique(roi);
i(find(i==0))=[];

fprintf('Finding centre points of ROIs for labels...');
for j = 1:length(i)
    vox    = find(roi==i(j));
    verts{j}  = v(vox,:);
    C(j,:) = spherefit(verts{j});
end
fprintf('  ... done! \n');
% % Plot the first roi, mark centre and label:
% scatter3(v(:,1),v(:,2),v(:,3),'k'); hold on
% scatter3(verts(:,1),verts(:,2),verts(:,3),'r')
% scatter3(C(:,1),C(:,2),C(:,3),'b*')

end

function [C,verts] = get_roi_centres0(pos,all_roi_tissueindex)
% Find centre points of rois
%
% starting at 0
v   = pos;
roi = all_roi_tissueindex;

i   = unique(roi);
%i(find(i==0))=[];
i = i(1):i(end); % don't skip empties
fprintf('Finding centre points of ROIs for labels...');
for j = 1:length(i)
    vox    = find(roi==i(j));
    verts{j}  = v(vox,:);
    C(j,:) = spherefit(verts{j});
end
fprintf('  ... done! \n');
% % Plot the first roi, mark centre and label:
% scatter3(v(:,1),v(:,2),v(:,3),'k'); hold on
% scatter3(verts(:,1),verts(:,2),verts(:,3),'r')
% scatter3(C(:,1),C(:,2),C(:,3),'b*')

end



function [Centre,A,B] = spherefit(X)
% Fit sphere to centre of vertices, return centre points
%
%

A =  [mean(X(:,1).*(X(:,1)-mean(X(:,1)))), ...
    2*mean(X(:,1).*(X(:,2)-mean(X(:,2)))), ...
    2*mean(X(:,1).*(X(:,3)-mean(X(:,3)))); ...
    0, ...
    mean(X(:,2).*(X(:,2)-mean(X(:,2)))), ...
    2*mean(X(:,2).*(X(:,3)-mean(X(:,3)))); ...
    0, ...
    0, ...
    mean(X(:,3).*(X(:,3)-mean(X(:,3))))];
A = A+A.';
B = [mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,1)-mean(X(:,1))));...
     mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,2)-mean(X(:,2))));...
     mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,3)-mean(X(:,3))))];
Centre=(A\B).';
end


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
















% Notes / Workings
%---------------------------------------------------
    %rotations - because x is orientated backward?
%     t  = 90;
%     Rx = [ 1       0       0      ;
%            0       cos(t) -sin(t) ;
%            0       sin(t)  cos(t) ];
%     Ry = [ cos(t)  0      sin(t)  ;
%            0       1      0       ;
%           -sin(t)  0      cos(t)  ];
%     Rz = [ cos(t) -sin(t) 0       ;
%            sin(t)  cos(t) 0       ;
%            0       0      1       ];
   %M = (Rx*(M'))';
   %M = (Ry*(M'))';
   %M = (Rz*(M'))';

   
