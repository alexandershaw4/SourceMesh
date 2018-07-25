function slice6()
% Take a figure with 1 axis and put it into a 6 subplots in a new figure.
% 
% For brains (see atemplate). 
%
% Rotates so three axes are L, topo, R, back, bottom & front
% 
% AS17

coloff;
f_c   = gcf;
f     = figure('Position', [622         195        1728        1139]);
% a(1)  = subplot(2,3,1,'align');
% a(2)  = subplot(2,3,2,'align'); 
% a(3)  = subplot(2,3,3,'align'); 
% a(4)  = subplot(2,3,4,'align');
% a(5)  = subplot(2,3,5,'align'); 
% a(6)  = subplot(2,3,6,'align'); 

a(1) = subplot('position',[0 2 1 1]/4); ... [left bottom width height]
a(2) = subplot('position',[1 2 1 1]/4);
a(3) = subplot('position',[2 2 1 1]/4);
a(4) = subplot('position',[0 1 1 1]/4);
a(5) = subplot('position',[1 1 1 1]/4);
a(6) = subplot('position',[2 1 1 1]/4);



axes_to_be_copied = findobj(f_c,'type','axes'); 
children = get(axes_to_be_copied,'children'); 
try
    overlay  = get(children{1}(end),'FaceVertexCData');
catch
    overlay  = get(children(end),'FaceVertexCData');
end

% set views:
v{1} = [270 0];  % L
v{2} = [0  90];  % Topo
v{3} = [90  0];  % R
v{4} = [0   0];  % back
v{5} = [0 -90];  % bottom
v{6} = [-180 0]; % front

% set aspects:
asp{1} = [2 3 2];
asp{2} = [2 3 2];
asp{3} = [2 3 2];
asp{4} = [2 3 2];
asp{5} = [2 3 2];
asp{6} = [2 3 2];

for k = 1:length(a)
    this = a(k);
    if iscell(children)
        for i = 1;
            stuff = handle2struct(children{i});    
            struct2handle(stuff,this);
            set(f,'currentaxes',this);
            set(gca,'visible','off');  
            axis(gca,'image');
        end
    else
            stuff = handle2struct(children);    
            struct2handle(stuff,this);
            set(f,'currentaxes',this);
            set(gca,'visible','off');   
            axis(gca,'image');
    end
    view(this,v{k});
    pbaspect(asp{k});
    axis(gca,'image');
end

for i = 1:length(a)
    p    = get(a(i), 'pos');
    p(3) = p(3) + 0.1;
    set(a(i), 'pos', p);
end
    