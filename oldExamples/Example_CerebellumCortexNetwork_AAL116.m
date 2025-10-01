
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
    'sourcemodel',{v vi}, 'method',{'user','spheres'},'optimise',0);


% Network - AAL116 (with cerebellum)
%-------------------------------------------------------
load AAL116 % load the AAL116 'source' positions 
N = zeros(116,116); % make an empty 116x116 network

% A cerebellar network (see 'labels' - nodes 91:108):
N = zeros(116,116);
N(91:108,91:108) = [ randi([0 1],18,18) .* randi([-6 6],18,18) ];


% just edges:
figure('position',[1000 249 885 729]);
atemplate('mesh','def3','network',N,'sourcemodel',{v vi})
view([0   17.5791])

% with nodes as well as edges:
figure('position',[1000 249 885 729]);
atemplate('mesh','def3','network',N,'sourcemodel',{v vi},'nodes',sum(N))
view([0   17.5791])
