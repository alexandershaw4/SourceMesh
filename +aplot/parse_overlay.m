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