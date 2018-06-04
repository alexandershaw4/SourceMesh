%Examples to use the InteractiveColorbar

%basic use
%
%Use the mouse pointer on the colorbar on the top and the bottom. When the
%mouse cursor changes in a hand, the CLim property of the image can be
%changed interactively. Use right mouse button on the InteractiveColorbar
%to set the Colormap.
figure;
imagesc(rand(100,100))
InteractiveColorbar;


%command line manipulation
figure
imagesc(rand(100,100))
ColorB=InteractiveColorbar;

%Set the CLim manually. Note setting the CLim property of the axes will not
%affect the colorbar. Set the CLim property on the InteractiveColorBar
%instead.
set(ColorB,'CLim',[0.25 0.75])

%Extend or limit the range of the InteractiveColorbar
set(ColorB,'CRange',[-1 2]); %Extend
set(ColorB,'CRange',[0.4 0.6]); %Limit

%Hide text labels
set(ColorB,'ShowCLimLabels',false)


%Advanced use: Connect multiple InteractiveColorbars
clear ColorB;figure;
for i=1:4    
     h(i)=subplot(2,2,i);
     imagesc(rand(100,100));
     axis off;
end

for i=1:4
    ColorB(i)=InteractiveColorbar(h(i));
end

%Link InterActiveColorbars:

%Move one Colorbar, all will follow
set(ColorB,'ConnectedColorBarObjects',ColorB);

%Disable connection for the first InteractiveColorbar
set(ColorB(1),'UpdateConnectedColorBars',false);




