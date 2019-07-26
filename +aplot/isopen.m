function data = isopen(data,in)
    % a wrapper on atemplate for generating the same plot for the left and
    % right hemis separately in different sublots
    s(1) = subplot(121);
    s(2) = subplot(122);
    
    data_in     = data;
    
    % Plot 1.
    in.hemi     = 'l';
    in.fighnd   = s(1);
    data        = sort_sourcemodel(data_in,in); % Sourcemodel vertices
    [mesh,data] = get_mesh(in,data);         % Get Surface
    [data,in]   = sort_template(data,in);    % Template space? (aal90/78/58)
    [mesh,data] = parse_mesh(mesh,in,data);  % Mesh to put stuff on
    data.mesh   = mesh;
    data        = parse_plots(data,in);      % Overlays, networks, etc
    data.in     = in;                        % Return everything
    data_l      = data;
    
    view([90 0]);
    
    % Plot 2.
    in.hemi     = 'r';
    in.fighnd   = s(2);
    data        = sort_sourcemodel(data_in,in); % Sourcemodel vertices
    [mesh,data] = get_mesh(in,data);         % Get Surface
    [data,in]   = sort_template(data,in);    % Template space? (aal90/78/58)
    [mesh,data] = parse_mesh(mesh,in,data);  % Mesh to put stuff on
    data.mesh   = mesh;
    data        = parse_plots(data,in);      % Overlays, networks, etc
    data.in     = in;                        % Return everything
    data_r      = data;
    
    view([270 0]);
    
    data   = [];
    data.l = data_l;
    data.r = data_r;

end