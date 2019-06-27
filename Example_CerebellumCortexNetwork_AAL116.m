
% Functional Overlay - AAL116 (with cerebellum)
%-------------------------------------------------------

% use atlas-reduced code: aal116
fun_data = zeros(116,1);
fun_data(91:108)=-8:9;   % activations in cerebellum only

% some cortical stuff too
fun_data(11) = -7;
fun_data(20) =  7;
fun_data(59) = -4;
fun_data(84) =  4;

% generate the figure this way:
figure('position',[1000 249 885 729]);
atemplate('mesh','def3','overlay',fun_data,'method',{'aal116','spheres'});

% or this way: supplying the sources and using explicit optimisation
load AAL116

figure('position',[1000 249 885 729]);
atemplate('mesh','def3','overlay',fun_data,...
    'sourcemodel',{v vi}, 'method',{'user','spheres'},'optimise',1);


% Network
%-------------------------------------------------------

load AAL116 % load the AAL116 'source' positions 

N = zeros(116,116); % make an empty 116x116 network

% 49  = L Sup Occ
% 105 = L Cerebellum 9

N(49,105) = 2;
N(105,49) = 2;

% just edges:
figure('position',[1000 249 885 729]);
atemplate('mesh','def3','network',N,'sourcemodel',{v vi})

% with nodes as well as edges:
figure('position',[1000 249 885 729]);
atemplate('mesh','def3','network',N,'sourcemodel',{v vi},'nodes',sum(N))