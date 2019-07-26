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
            [mesh, data] = aplot.convert_mesh(mesh,data);
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
                v(:,i) = aplot.rescale(v(:,i),B(:,i));
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
            if strcmp(lower(mesh),'def5')
                nmesh         = load('ft_mesh_mni');
                mesh          = [];
                mesh.vertices = nmesh.vertices;
                mesh.faces    = nmesh.faces;
            end
    end
    try
        close(wb);  
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
                v(:,i) = aplot.rescale(v(:,i),B(:,i));
            end
            
            % return scaled mesh
            mesh            = [];
            mesh.faces      = fv.faces;
            mesh.vertices   = v;
            data.mesh       = mesh;    
            
            close(wb);        
end

end