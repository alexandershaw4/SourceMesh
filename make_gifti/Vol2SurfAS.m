function g = Vol2SurfAS(file,type,varargin)
% Convert ctf mri or (fieldtrip) matlab format mri to a beautiful, smooth 
% isosurface gifti object. Designed for use with atemplate plot tool.
%
% usage: g = Vol2SurfAS(file,type)
%        g = Vol2SurfAS('myctfmri.mri','ctf')
%        g = Vol2SurfAS('mysegmri.mat','ft')
%
%
% optional paired additional inputs:
% 'smooth', n
%
%
% AS

smth   = 0;
for i  = 1:length(varargin)
    if strcmp(varargin{i},'smooth'); smth = varargin{i+1}; end
end

switch type
    case {'ctfmri','ctf'}
        mri = ft_read_mri('010213-50_coreg.mri','dataformat', 'ctf_mri4');
    case {'ft','matlab'}
        load(file);
        %file = 'coregseg_0001.mat';
end


% normalise & isosurface 
cfg.template   = [fileparts(mfilename('fullpath')) '/T1.nii'];
cfg.spmversion = 'spm12';
mri            = ft_volumenormalise(cfg,mri);
mri            = ft_volumereslice([], mri);
segmentedmri   = ft_volumesegment(cfg, mri);
V              = isosurface(segmentedmri.gray,.5);

% smooth and centre
dV = sms(V.vertices,V.faces,20,.5);
cV = dV - repmat(spherefit(dV),[size(dV,1),1]);

v.faces    = V.faces;
v.vertices = cV;

t.vertices = [ v.vertices(:,2), v.vertices(:,1), v.vertices(:,3) ];
t.faces    = [ v.faces(:,2),    v.faces(:,1),    v.faces(:,3) ];

g = gifti(t);

if smth > 0
    [N.faces, N.vertices] = reducepatch(g.faces,g.vertices,smth);
    g = gifti(N);
end


% % smooth [less] and centre // quicker:
% smV = sms(V.vertices,V.faces,1,4);
% smV = smV - repmat(spherefit(smV),[size(smV,1),1]);
% smth.faces    = V.faces;
% smth.vertices = smV;
% sm            = gifti(smth);



















% % rotate about z
% Q     = @squeeze;
% for i = 1:size(v.vertices,1)
%     newv(i,:) = Rz*Q(v.vertices(i,:))';   
% end
% 
% new.vertices = newv;
% new.faces = v.faces;
% newg = gifti(new);
% plot(newg)

% Notes / Workings
%---------------------------------------------------
    %rotations - because x is orientated backward?
%     t  = 315;
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






% [no,el,regions,holes] = v2s(segmentedmri.gray,0.2,5);
% 
% c  = spherefit(no);
% no = no - repmat(c,[size(no,1),1]);
% 
% e.vertices = no;
% e.faces    = el;
% V          = gifti(e);

% % smooth surface
% newnode=sms(no,el,20,.5);
% 
% sm.vertices = newnode;
% sm.faces = el;
% smooth = gifti(sm);


% % surface to mesh?
% [node,elem,face]=s2m(no,el,1,1);
% 
% smesh.vertices = node;
% smesh.faces = face;
% gmesh = gifti(smesh);