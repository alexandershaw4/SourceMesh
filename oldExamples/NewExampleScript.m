% An example script demonstrating how to use atemplate plot function
%
% Generates a random AAL90 network and places in left subplot
% Generates a functional overlay vector and plots in right subplot
% Links the 2 images, so that they rotate together


% open a big figure window 
%--------------------------------------------------------------------------
figure('position',[1000 142 1343 836]);


% make a spoof AAL90 network and plot on the left
%--------------------------------------------------------------------------
net                                        = zeros(90,90);
net(randi([1 90],20,1),randi([1 90],20,1)) = randi([-7 7],20,20);

s(1) = subplot(121); atemplate('network',net,'fighnd',s(1));


% make AAL90 'functional overlay' and plot on the right
%--------------------------------------------------------------------------
over = randi([-7 9],90,1);

s(2) = subplot(122); atemplate('overlay',over(:),'fighnd',s(2),'nocolbar');


% use linksubplots function to link rotation, so the two images move
% together in the window
%--------------------------------------------------------------------------
linksubplot(s)



% Alternatively, plot overlay in left and right subplots, linked:
%--------------------------------------------------------------------------
figure('position',[1000 142 1343 836]);

s(1) = subplot(121); atemplate('overlay',over,'hemi','l','fighnd',s(1),'nocolbar')
s(2) = subplot(122); atemplate('overlay',over,'hemi','r','fighnd',s(2),'nocolbar')

linksubplots(s)