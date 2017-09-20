function slice3()
% Take a figure with 1 axis and put it into a 3 subplots in a new figure.
% 
% For brains (see atemplate). Rotates so three axes are L, topo & R.
% 
% AS17

coloff;
f_c   = gcf;
f     = figure;
a(1)  = subplot(1,3,1);
a(2)  = subplot(1,3,2); 
a(3)  = subplot(1,3,3); 

axes_to_be_copied = findobj(f_c,'type','axes'); 
children = get(axes_to_be_copied,'children'); 
try
    overlay  = get(children{1}(end),'FaceVertexCData');
catch
    overlay  = get(children(end),'FaceVertexCData');
end

% set views:
v{1} = [270 0]; % L
v{2} = [0  90]; % Topo
v{3} = [90  0]; % R

asp{1} = [2 3 2];
asp{2} = [1 1 1];
asp{3} = [2 3 2];

for k = 1:3
    this = a(k);
    if iscell(children)
        for i = 1;
            stuff = handle2struct(children{i});    
            struct2handle(stuff,this);
            set(f,'currentaxes',this);
            set(gca,'visible','off');        
        end
    else
            stuff = handle2struct(children);    
            struct2handle(stuff,this);
            set(f,'currentaxes',this);
            set(gca,'visible','off');        
    end
    view(this,v{k});
    pbaspect(asp{k});
end

set(f, 'Position', [100, 100, 2300, 600]);

for i = 1:3
    p    = get(a(i), 'pos');
    p(3) = p(3) + 0.1;
    set(a(i), 'pos', p);
end
    

