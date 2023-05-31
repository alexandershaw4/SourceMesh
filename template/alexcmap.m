function map = alexcmap(m)

[map2,map3,map4]=cubric_meg_palettes;

map(:,1) = spm_vec([map3(:,1) map3(:,1)]');
map(:,2) = spm_vec([map3(:,2) map3(:,2)]');
map(:,3) = spm_vec([map3(:,3) map3(:,3)]');

if nargin > 0

    for i = 1:3
        %im(:,i) = interp1((1:256)-(256/2),map(:,i),(1:m)-(m/2),'cubic');
        im(:,i) = resample(map(:,i),m,256);
    end
    map = im;
    map = max(min(round(map*100)/100,1),0);

end
