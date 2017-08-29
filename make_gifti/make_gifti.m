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


properties 
    mri
    type
    g
    i
end

methods
    
    function obj = make_gifti(mri, type)
             obj.mri  = mri;
             obj.type = type;
    end
    function obj = makesurf(obj, smth)
    if nargin < 2; smth = .15; end
        obj.g = Vol2SurfAS(obj.mri,obj.type,smth);
    end
    
    function obj = reduce(obj,smth)
    if nargin < 2; smth = .15; end
        [obj.g.faces, obj.g.vertices] = reducepatch(obj.g.faces,obj.g.vertices,smth);
    end
    
    function obj = inflate(obj,howmuch)
    if nargin < 2; howmuch = 400; end
            obj.i = spm_mesh_inflate(obj.g,howmuch);
    end
    function obj = smooth(obj,smth)
    if nargin < 2; smth = 0.5; end
        dV = sms(V.vertices,V.faces,5,smth);
        obj.g.vertices = dV;
    end
    
    
end


end