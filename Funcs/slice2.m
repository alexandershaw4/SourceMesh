function slice2()
% Take a figure with 1 axis and put it into a 2 subplots in a new figure.
% 
% For brains (see atemplate). Rotates so three axes are L, R.
% 
% AS17

coloff;
f_c   = gcf;
f     = figure;
a(1)  = subplot(1,2,1);
a(2)  = subplot(1,2,2);  

axes_to_be_copied = findobj(f_c,'type','axes'); 
children = get(axes_to_be_copied,'children'); 
overlay  = get(children{1}(end),'FaceVertexCData');

% set views:
v{1} = [270 0]; % L
v{2} = [90  0]; % R

asp{1} = [3 5 3];
asp{2} = [3 5 3];

for k = 1:2
    this = a(k);
    for i = 1;
        stuff = handle2struct(children{i});    
        struct2handle(stuff,this);
        set(f,'currentaxes',this);
        camlight('headlight')
        set(gca,'visible','off');        
    end
    view(this,v{k});
    pbaspect(asp{k});
end

set(f, 'Position', [100, 100, 1800, 600]);

for i = 1:2
    p    = get(a(i), 'pos');
    p(3) = p(3) + 0.1;
    set(a(i), 'pos', p);
end


