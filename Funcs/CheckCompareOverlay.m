function t = CheckCompareOverlay(D)
% give the output structure of an atemplate call containing 'overlay'
% structure

Name   = {'Orig Data' 'Overlay Model' 'Sq Error'}';
Mean   = [mean(D.overlay.orig) mean(D.overlay.data)]';
Mode   = [mode(D.overlay.orig) mode(D.overlay.data)]';
Median = [median(D.overlay.orig) median(D.overlay.data)]';
Min    = [min(D.overlay.orig) min(D.overlay.data)]';
Max    = [max(D.overlay.orig) max(D.overlay.data)]';

Mean(3) = (Mean(1)-Mean(2)).^2;
Mode(3) = (Mode(1)-Mode(2)).^2;
Median(3) = (Median(1)-Median(2)).^2;
Min(3)  = (Min(1)-Min(2)).^2;
Max(3)  = (Max(1)-Max(2)).^2;

t = table(Name,Mean,Mode,Median,Min,Max);