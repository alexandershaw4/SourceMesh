classdef make_gifti < handle
% Constructor for make_gifti / Vol2SurfAS
%
% Convert ctf mri (.mri) to cortical surface
%
% use:
% this = make_gifti('my-mri.mri','ctf'); % initiate object
% this.makesurf;                         % run segment and surface function
% this.inflate;                          % inflate the resulting mesh
% this.reduce;                           % simplify mesh
% this.smooth;                           % smooth mesh
% 
% use:
% this = make_gifti([],'spm'); % initiate object
% this.makesurf;               % load spm 8k vertex surface
% this.smooth / inflate 
%
% AS

properties 
    mri
    type
    g
    i
    orig
end

methods
    
    function obj = make_gifti(mri, type)
             obj.mri  = mri;
             obj.type = type;
    end
    function obj = makesurf(obj, smth)
    if nargin < 2; smth = .15; end
        switch obj.type
            case 'ctf';
                obj.g    = Vol2SurfAS(obj.mri,obj.type,smth);
                obj.orig = obj.g;
            case {'gii','gifti'}
                obj.g    = gifti(obj.mri);
                obj.orig = obj.g;
            case 'spm'
                obj.g    = gifti('spm.surf.gii');
                obj.orig = obj.g;
        end
    end
    
    function obj = reduce(obj,smth)
    if nargin < 2; smth = .15; end
        [obj.g.faces, obj.g.vertices] = reducepatch(obj.g.faces,obj.g.vertices,smth);
    end
    
    function obj = inflate(obj,howmuch)
    if nargin < 2; 
        switch obj.type;
            case 'ctf'; howmuch = 400; 
            case 'spm'; howmuch = 60;
        end
            obj.i = spm_mesh_inflate(obj.g,howmuch);
    end
    end
    function obj = smooth(obj,smth)
    if nargin < 2; smth = 0.5; end
        dV = sms(obj.g.vertices,obj.g.faces,5,smth);
        obj.g.vertices = dV;
    end
    function obj = recentre(obj)
        V  = obj.g.vertices;
        cV = V - repmat(spherefit(V),[size(V,1),1]); 
        obj.g.vertices = cV;
    end
    
    
end
end