function varargout = plot(varargin)
% plot method for GIfTI objects
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: plot.m 5888 2014-02-19 19:54:12Z guillaume $

% if ishandle(varargin{1})
%     h = figure(varargin{1});
% else
%     h = figure;
%     %axis equal;
%     %axis off;
%     %camlight;
%     %camlight(-80,-10);
%     %lighting phong;
% end
% cameramenu;


cdata = [];
ax = [];
indc = 0;
if nargin == 1
    this = varargin{1};
    h = gcf;
else
    if ishandle(varargin{1})
        ax = varargin{1};
        h = figure(get(ax,'parent'));
        this = varargin{2};
        
        % allow pass of cdata as varargin{3}
        if nargin == 3;
            cdata = subsref(varargin{3},struct('type','.','subs','cdata'));
            indc  = 1;
        end
    elseif strcmp(varargin{2},'fighnd')
        ax = varargin{3};
        h = figure(get(ax,'parent'));
        this = varargin{1};
        
        % allow pass of cdata as varargin{3}
        if nargin == 4;
            cdata = subsref(varargin{4},struct('type','.','subs','cdata'));
            indc  = 1;
        end        
    else
        this = varargin{1};
        h = gcf;
        cdata = subsref(varargin{2},struct('type','.','subs','cdata'));
    end
    % changed in light of above
    if nargin > 2 && isempty(indc)
        indc = varargin{3};
    else
        indc = 1;
    end
end

if isempty(ax), ax = axes('Parent',h); end

% Alex add: only make size x=y if not a subplot (i.e. no handle)
if ~ishandle(varargin{1});
    axis(ax,'equal');
end

axis(ax,'off');
hp = patch(struct(...
    'vertices',  subsref(this,struct('type','.','subs','vertices')),...
    'faces',     subsref(this,struct('type','.','subs','faces'))),...
    'FaceColor', 'b',...
    'EdgeColor', 'none',...
    'Parent',ax);

if ~isempty(cdata)
    set(hp,'FaceVertexCData',cdata(:,indc), 'FaceColor','interp')
end

axes(ax);
camlight;
camlight(-80,-10);
lighting phong;
axes(ax);
cameramenu;

if nargout
    varargout{1} = hp;
end
