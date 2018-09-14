function varargout = CheckOrientationMesh(varargin)
% CHECKORIENTATIONMESH MATLAB code for CheckOrientationMesh.fig
%      CHECKORIENTATIONMESH, by itself, creates a new CHECKORIENTATIONMESH or raises the existing
%      singleton*.
%
%      H = CHECKORIENTATIONMESH returns the handle to a new CHECKORIENTATIONMESH or the handle to
%      the existing singleton*.
%
%      CHECKORIENTATIONMESH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHECKORIENTATIONMESH.M with the given input arguments.
%
%      CHECKORIENTATIONMESH('Property','Value',...) creates a new CHECKORIENTATIONMESH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CheckOrientationMesh_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CheckOrientationMesh_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CheckOrientationMesh

% Last Modified by GUIDE v2.5 14-Sep-2018 15:24:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CheckOrientationMesh_OpeningFcn, ...
                   'gui_OutputFcn',  @CheckOrientationMesh_OutputFcn, ...
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


% --- Executes just before CheckOrientationMesh is made visible.
function CheckOrientationMesh_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CheckOrientationMesh (see VARARGIN)

% Choose default command line output for CheckOrientationMesh
handles.output = hObject;

S = varargin{1};
handles.input = inputname(1);
mesh.vertices = double(S.vertices);
mesh.faces    = double(S.faces);

axes(handles.axes1);
p = patch('faces',mesh.faces,'vertices',mesh.vertices);

C    = [.5 .5 .5];

set(p,'FaceColor',C);
box off;
grid off; 
set(p,'EdgeColor','none')
set(gca,'visible','off');
lighting phong
camlight left
camlight right

handles.mesh = mesh;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CheckOrientationMesh wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CheckOrientationMesh_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
varargout{1} = handles.mesh;

delete(handles.figure1);

% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1

% flip top / bottom
comps = allchild(handles.axes1);
p = comps(3);
v = p.Vertices;
v(:,3) = v(:,3)*-1;
set(p,'Vertices',v);

set(hObject,'Value',0)

% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2

% flip L/R
comps = allchild(handles.axes1);
p = comps(3);
v = p.Vertices;
v(:,1) = v(:,1)*-1;
set(p,'Vertices',v);

set(hObject,'Value',0)

% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3

% flip front / back
comps = allchild(handles.axes1);
p = comps(3);
v = p.Vertices;
v(:,2) = v(:,2)*-1;
set(p,'Vertices',v);

set(hObject,'Value',0)

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Export / Done
comps = allchild(handles.axes1);
p = comps(3);
v = p.Vertices;
f = p.Faces;

mesh.vertices = v;
mesh.faces = f;
%assignin('base','mesh',mesh);
%assignin('caller','mesh',mesh);
handles.mesh = mesh;
guidata(hObject, handles);
%CheckOrientationMesh_OutputFcn(hObject, eventdata, handles) 



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
%delete(hObject);

if isequal(get(hObject, 'waitstatus'), 'waiting')
% The GUI is still in UIWAIT, us UIRESUME
uiresume(hObject);
else
% The GUI is no longer waiting, just close it
delete(hObject);
end