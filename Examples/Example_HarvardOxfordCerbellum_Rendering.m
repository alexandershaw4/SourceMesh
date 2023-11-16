
% Render Harvard-Oxford data with (AAL) Cerebellum 
%-------------------------------------------------------

% Has 74 nodes: 48 HarvOx Atlas, 26 AAL Cerebellar/Vermis
load HOA_Cerebellum.mat

fun_data = zeros(74,1);

% see 'labels':
% vals 1:48  = HOA regions
%      49:74 = AAL cerebellar
%

fun_data(1:48)  = randi([0 7],48,1); 
fun_data(49:74) = -7;

% generate overlay figure this way:
figure('position',[1000 249 885 729]);
atemplate('mesh','def3','overlay',fun_data,'method',{'hoac','spheres'});

% or this way: supplying the sources and using explicit optimisation
load HOA_Cerebellum.mat

figure('position',[1000 249 885 729]);
atemplate('mesh','def3','overlay',fun_data,...
    'sourcemodel',{v vi}, 'method',{'user','spheres'});

