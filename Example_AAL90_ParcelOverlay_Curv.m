% Example script: generate cortical mesh of AAL90 parcellation with
% curature and inflation

% Make reduced version of default mesh 5
m  = atemplate('mesh','def5');close;drawnow;
m  = maker(m.mesh.faces,m.mesh.vertices);
fv = reducepatch(m,.2);

% AAL90 'Overlay' data:
z = zeros(1,90);

z([19 20]) = [-5 5];   % motor
z([43 44]) = [-10 10]; % calcarine

 afigure; 
 D = atemplate('inflate',...
        'overlay',{'curvature',z},...
        'method' ,{'aal_light','spheres'} ,...
        'thresh' ,0, ...
        'mesh'   ,fv );
    