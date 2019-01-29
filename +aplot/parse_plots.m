function data = parse_plots(data,i)

% unpack triggers
inputs = i;

% overlays
if isfield(inputs,'L')
    % copy over overlay options
    data.overlay.peaks      = i.peaks;
    data.overlay.components = i.components;
    data.overlay.pca        = i.pca;
    data.overlay.method     = i.method;
    data.overlay.depth      = i.depth;
    data.overlay.tf_interactive = i.tf_interactive;
    
    if isfield(i,'funcaffine')
        data.overlay.affine = i.funcaffine;
    end
    
    data = overlay(data, (i.L),i.write,i.fname,i.colbar);
end 

isover = exist('L','var') || exist('V','var');
if  isover && exist('A','var') 
    i.colbar = 0;
    alpha(.2);
end

% networks
if isfield(inputs,'A')
    data = connections(data,i.A,i.colbar,i.write,i.fname,i.netcmap); 
end 

% tracts
if isfield(inputs,'T')
    data = drawtracks(data,i.T,i.H);                  
end 

% nodes
if isfield(inputs,'N')
    data = drawnodes(data, i.N);                 
end 

% labels
data = parse_labels(i,data);

% video
if isfield(inputs,'V')
    tv = 1:size(i.V,2);
    try tv = i.times; end
    data = video(data,i.V,1,i.fpath,tv); 
end


end