function varargout = atemplate_gui(varargin)
% For plotting networks, overlays, nodes, labels & tracks* derived from the
% AAL90 atlas.
%
% Use the menus to load the things you need and click update plot.
% 
% This is a simple frontend for my function atemplate.m
% 
% 
% 
% 
% AS
% See also: atemplate

% Edit the above text to modify the response to help atemplate_gui

% Last Modified by GUIDE v2.5 10-Aug-2017 14:54:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @atemplate_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @atemplate_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before atemplate_gui is made visible.
function atemplate_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to atemplate_gui (see VARARGIN)

% Choose default command line output for atemplate_gui
handles.output = hObject;

% initialise everything
handles.network = 0;
handles.overlay = 0;
handles.nodes   = 0;
handles.labels  = 0;
handles.tracks  = 0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes atemplate_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = atemplate_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function png_Callback(hObject, eventdata, handles)
% hObject    handle to png (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName,FilterIndex] = uiputfile({'*.png'},'Filename');

print(gcf,[PathName FileName],'-dpng','-r600');

% --------------------------------------------------------------------
function fig_Callback(hObject, eventdata, handles)
% hObject    handle to fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName,FilterIndex] = uiputfile({'*.fig'},'Filename');

savefig([PathName FileName]);

% --------------------------------------------------------------------
function NodeEdge_Callback(hObject, eventdata, handles)
% hObject    handle to NodeEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName,FilterIndex] = uigetfile({'*.edge'},'Select edge or node file');
A = dlmread([PathName '/' FileName]);

handles.network = 1;
handles.edges   = A;

L = menu('Labels?','Yes','No');
% popup = uicontrol('Style', 'popup',...
%            'String', {'Labels','No Labels'},...
%            'Position', [20 340 100 50]);

if L == 1;
    handles.labels = 1;
else
    handles.labels = 0;
end

guidata(hObject, handles);


% --------------------------------------------------------------------
function Overlay_Callback(hObject, eventdata, handles)
% hObject    handle to Overlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
W = menu('Load from where','Matfile','Matlab workspace');

if W == 2;
    list = evalin('base','whos');
    opt  = menu('Select variable (1x90)',{list.name});
    list = {list.name};
    var  = list{opt};
    handles.O = evalin('base',var);
    handles.overlay = 1;
end
if W == 1;
    [FileName,PathName,FilterIndex] = uigetfile({'*.mat'},'Select matfile with variable O');
    x = load([PathName FileName]);
    handles.O = x.O;
    handles.overlay = 1;
end
guidata(hObject, handles);
    

% --------------------------------------------------------------------
function Nodes_Callback(hObject, eventdata, handles)
% hObject    handle to Nodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Tracks_Callback(hObject, eventdata, handles)
% hObject    handle to Tracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

atemplate;

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% update the plot - i.e. call atemplate
str = {};  cla;

if handles.network;
    str = {str, 'network', handles.edges};
end
if handles.labels;
    str = {str{:}, 'labels'};
end
if handles.nodes;
    str = {str{:}, 'nodes', handles.N};
end
if handles.overlay;
    str = {str{:}, 'overlay',handles.O};
end
if handles.tracks;
    str = {str{:}, 'tracks',handles.T,handles.H};
end
if isempty(str{1});
    str = str(2:end);
end
atemplate(str{:});
