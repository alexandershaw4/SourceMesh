classdef InteractiveColorbar < hgsetget
    properties (Constant)
        %Default properties for UI
        COLORMAP_IMAGE_WIDTH=5
        COLORBAR_HEIGHT_RATIO=0.95;
        COLORBAR_WIDTH=20;
        DEFAULT_COLORMAP='jet';
        RULER_COLOR='k'
        RULER_POINT_SHAPE='none'
        LINE_WIDTH=3;
        Type='InteractiveColorbar'
        NUMTEXT_FORMAT='%10.1e'
        
    end
    
    properties (Access=private) 
        %Handles to components
        AxesHandle
        ColorBarHandle
        LabelHandle
        PointHandle
        LineHandle
        
        %The matlab colomaps of 6bit (64 steps) are interpolated to 1024
        %steps (10 bit)
        ColorMapSteps=1024; %Steps for displaying the colorbar image
        OffsetFrac
        %Stored Values for dependend properties
        StoredCRange
        StoredCLim
        StoredConnectedColorBarObjects
        
        %Listeners for updates in connected colorbars (see
        %ConnectedColorBarObjects)
        ColorBarListeners
        
        %Listener for resizes of axes, colorbar needs resizing then
        ResizeListener
        
        %Listen for deletion of axes, delete colorbar if so        
        DeleteListener
        
        %Colormaps for context menu (right click on colorbar)
        ColorMaps={'Gray','Inverted Gray','Jet','HSV','Hot','Cool','Spring','Summer','Autumn','Winter','Bone','Copper','Pink','Lines','Colormap Editor'}
        
        %Handle for the context menu
        ContextMenu
    end
    
    properties(Dependent,Access=public)
        %Cell array of colorBar objects, changes in this colorbar object
        %will update the connected Colorbars
        ConnectedColorBarObjects
        
        %Show or Hide Value Labels next to the colorbar
        ShowCLimLabels
        
        %Set the colormap
        ColorMap
        
        %Set the CLim, Update is triggered and CLim for the connected axes
        %is set
        CLim
                
        %Set the valid range for CLim Values, Update is triggered and
        %ColorBar may stretch/shrink to enable the new Range
        CRange
        
    end
    
    properties(Dependent,SetAccess=public,GetAccess=private)
        CLimWithoutUpdate
        CLimWithoutUpdateConnectedColorbarObjects
        CRangeWithoutUpdateConnectedColorbarObjects
        CRangeWithoutUpdate        
    end
    
    properties (Dependent,Access=private)
       
       AlphaDataForColorBarImage %Needed for the UI
       ColorBarImage %Image (visible gradient) on the colorbar
       ColorBarPosition %4 vector for position calculated based on the obj.AxesHandle position property
       MouseCLim %Read Mouse position and calulate new CLim value
    end
    
    properties(Access=public)
        UpdateConnectedColorBars=false; %Trigger an update for connect colorbars
    end
        
    events (ListenAccess=public, NotifyAccess=private)
        UpdateByUser %Triggered when user interacts with the colorbar, external objects may listen for this event
    end
        
    methods 
        
        %Constructor
        function obj=InteractiveColorbar(AxesHandle)
            
            %apply on current axes if no axes is specified
            if nargin<1;AxesHandle=gca;end
            
            %Check if handle is an AxesHandle
            assert(ishandle(AxesHandle),'No valid handle as input');
            assert(strcmp(get(AxesHandle,'Type'),'axes'),'Handle is not an axes handle');
                                                                    
            %Some fail safe methods to prevent adding a colorbar to a
            %colorbar. Or adding a colorbar when a colorbar is already
            %present. Information is stored using set appdata at the end of
            %the constructor.
            axestype=getappdata(AxesHandle,'ColorBarType');            
            if isempty(axestype);axestype='empty';end
            
            ShrinkAxesWidth=2.*obj.COLORBAR_WIDTH;
            
            switch axestype
                case 'InteractiveColorbar'
                    newAxesHandle=getappdata(AxesHandle,'AxesHandle');
                    delete(AxesHandle);
                    AxesHandle=newAxesHandle;  
                    ShrinkAxesWidth=0;
                case 'AxesWithInteractiveColorbar'
                    ColorBHandle=getappdata(AxesHandle,'ColorBarHandle');
                    if ishandle(ColorBHandle);delete(ColorBHandle);end
                    ShrinkAxesWidth=0;
            end
            
            %Store axes handle
            obj.AxesHandle=AxesHandle;
            set(obj.AxesHandle,'Units','pixels');
            AxesPos=get(obj.AxesHandle,'Position');           
            AxesCLim=get(AxesHandle,'CLim');
            obj.StoredCLim=AxesCLim;
            obj.StoredCRange=AxesCLim;                                                           
            
            %Determine position to place the colorbar and resize the axes
            
            AxesPos(3)=AxesPos(3)-ShrinkAxesWidth;
            if all(AxesPos>=0)
                set(obj.AxesHandle,'Position',AxesPos);
            else
                warning('ColorBar may overlap axes');
            end
                                               
            obj.CreateColorBar
            
            %Set pointer behaviour such that the pointer is a hand on the
            %movable points
            iptPointerManager(gcf);
            iptSetPointerBehavior(obj.PointHandle,@(obj,src,event) set(gcf,'Pointer','hand'));
            iptSetPointerBehavior(obj.LineHandle,@(obj,src,event) set(gcf,'Pointer','hand'));
                                                           
            %This sets everything to the correct position
            obj.Update;
            
            %Set default colormap
            obj.ColorMap=obj.DEFAULT_COLORMAP;                                    
            obj.ResizeListener=addlistener(obj.AxesHandle,'Position','PostSet',@obj.Resize);
            set(obj.ColorBarHandle,'DeleteFcn',@(src,event) delete(obj));
        end
     
        
        %Destructor
        function delete(obj)
            obj.ConnectedColorBarObjects=[];                        
            delete(obj.ColorBarListeners);
            delete(obj.ResizeListener);
            delete(obj.DeleteListener);
        end
                
        %Setters and getters        
        function Range=get.CRange(obj)
            Range=obj.StoredCRange;
        end
        
        function set.CRange(obj,Range)           
            obj.SetCRange(Range,true);
        end
        
        function set.CRangeWithoutUpdate(obj,Range)
            obj.SetCRange(Range,false);
        end
        
         function set.CRangeWithoutUpdateConnectedColorbarObjects(obj,Range)
            obj.SetCRange(Range,false)
            obj.Update(false);            
        end
        
        
        function SetCRange(obj,Range,Update)
             assert(isnumeric(Range) && numel(Range)==2 && Range(2)>Range(1),'Range should be a two element vector and element 2 should be greater than element 1');                                             
           assert(~all(Range<obj.CLim(1)) && ~all(Range>obj.CLim(2)),'Range must have overlap with CLim');
           if all(Range==obj.StoredCRange);return;end
                                                                  
           newCLim=obj.CLim;
           if Range(1)>obj.CLim(1)
               newCLim(1)=Range(1);
           end
           if Range(2)<obj.CLim(2)
               newCLim(2)=Range(2);
           end
           
           obj.StoredCRange=Range;
           obj.CLimWithoutUpdate=newCLim;
           if Update
            obj.Update;
           end
        end
            
            
        
        function ColorBarImage=get.ColorBarImage(obj)
           %Create colormap image which is displayed on the bar
            StepSize=1/(obj.ColorMapSteps);            
            ColorBarImage=(0:StepSize:1)';
            ColorBarImage=repmat(ColorBarImage,1,obj.COLORMAP_IMAGE_WIDTH); 
                                   
        end
        
        function AlphaData=get.AlphaDataForColorBarImage(obj)
            %Alpha data to make left and right transparent, here the
            %movable points are plotted
            AlphaData=ones(size(obj.ColorBarImage));
            AlphaData(:,1)=0;
            AlphaData(:,end)=0;
           
        end
        
        function Position=get.ColorBarPosition(obj)
            AxesPos=get(obj.AxesHandle,'Position');
            
            Position(1)=AxesPos(1)+AxesPos(3)+0.5.*obj.COLORBAR_WIDTH;
            Position(2)=AxesPos(2)+0.5*(1-obj.COLORBAR_HEIGHT_RATIO).*AxesPos(4);
            Position(3)=obj.COLORBAR_WIDTH;
            Position(4)=AxesPos(4).*obj.COLORBAR_HEIGHT_RATIO;
        end
        
        function set.ShowCLimLabels(obj,Value)
            if Value
                set(obj.LabelHandle,'Visible','on');
            else
               set(obj.LabelHandle,'Visible','off');
            end
        end
        
        function set.ConnectedColorBarObjects(obj,ColorBarObjects)
            if ~isempty(ColorBarObjects)
                assert(isa(ColorBarObjects,'InteractiveColorbar'),'The right hand side of the assignment should be a InteractiveColorbarObject or a cell array of InteractiveColorbarObjects');                
            end
            
            %Remove self 
            ColorBarObjects(eq(ColorBarObjects,obj))=[];
            
            obj.StoredConnectedColorBarObjects=ColorBarObjects;
            if ~isempty(obj.StoredConnectedColorBarObjects)
                obj.UpdateConnectedColorBars=true;
            else
                obj.UpdateConnectedColorBars=false;
            end
        end
        
        function ColorBarObjects=get.ConnectedColorBarObjects(obj)
            ColorBarObjects=obj.StoredConnectedColorBarObjects;
            %Remove possible deleted objects
            if ~isempty(ColorBarObjects)
                ColorBarObjects=ColorBarObjects(isvalid(ColorBarObjects));
            end
        end
                
        function set.ColorMap(~,ColorMapStr)
            switch lower(ColorMapStr)
                case 'inverted gray'
                    CMap=1-(0:1/63:1)';
                    CMap=repmat(CMap,1,3);
                    colormap(CMap)                
                otherwise
                    colormap(ColorMapStr);
            end
        end
        
        function Lim=get.MouseCLim(obj)
            pos=get(obj.ColorBarHandle,'CurrentPoint');
            Range=obj.CRange;
                       
            LimFrac=pos(1,2)./(obj.ColorMapSteps);           
            Lim=LimFrac.*(Range(2)-Range(1))+Range(1);        
        end
                              
        function set.CLim(obj,Lim)
          obj.SetCLim(Lim,true);
          
        end
        
        function set.CLimWithoutUpdate(obj,Lim)
            obj.SetCLim(Lim,false);          
        end
        
        function set.CLimWithoutUpdateConnectedColorbarObjects(obj,Lim)
            obj.SetCLim(Lim,false)
            obj.Update(false);            
        end
        
            
        
        function SetCLim(obj,Lim,Update)
              assert(isnumeric(Lim) && numel(Lim)==2,'Set value for CLim is invalid') 
            if all(Lim==obj.StoredCLim);return;end
            if Lim(1)>=Lim(2);return;end
            
            Range=obj.CRange;
            if Lim(1) < Range(1);Lim(1)=Range(1);end
            if Lim(2) > Range(2);Lim(2)=Range(2);end
            
            obj.StoredCLim=Lim;
            if Update
                obj.Update;
            end
        end
        
        function Lim=get.CLim(obj)
            Lim=obj.StoredCLim;
        end
                
        
    end
    
    methods (Access=public)                                    
        function Resize(obj,~,~)
            if ~ishandle(obj.ColorBarHandle);return;end
            set(obj.ColorBarHandle,'Position',obj.ColorBarPosition);
        end                
    end
    
    
    methods (Access=private)
        %Create UI Elements
        
        function CreateColorBar(obj)
            axes(obj.AxesHandle);
            obj.ColorBarHandle=axes('Units','pixels','Position',obj.ColorBarPosition);
            set(obj.ColorBarHandle,'Color',get(gcf,'Color'));
            
            %Put colormap image on the axes
            Imageh=imagesc(obj.ColorBarImage);
            set(Imageh,'AlphaData',obj.AlphaDataForColorBarImage);
            axis xy;axis off;
            set(obj.ColorBarHandle,'Clim',obj.CLim);
            
            %Create movable points at the sides
            obj.PointHandle(1)=line(0.5,obj.ColorMapSteps.*obj.CLim(1)+0.5);
            obj.PointHandle(2)=line(obj.COLORMAP_IMAGE_WIDTH+0.5,obj.ColorMapSteps.*obj.CLim(1)+0.5);
            
            obj.PointHandle(3)=line(0.5,obj.ColorMapSteps.*obj.CLim(2)+0.5);
            obj.PointHandle(4)=line(obj.COLORMAP_IMAGE_WIDTH+0.5,obj.ColorMapSteps.*obj.CLim(2)+0.5);
            
            obj.LineHandle(1)=line([0.5 obj.COLORMAP_IMAGE_WIDTH+0.5],[1 1].*obj.ColorMapSteps.*obj.CLim(1)+0.5);
            obj.LineHandle(2)=line([0.5 obj.COLORMAP_IMAGE_WIDTH+0.5],[1 1].*obj.ColorMapSteps.*obj.CLim(1)+0.5);
            
            
            
            set(obj.PointHandle,'Marker',obj.RULER_POINT_SHAPE)
            set(obj.PointHandle,'MarkerFaceColor',obj.RULER_COLOR,'Color',obj.RULER_COLOR);
            set(obj.PointHandle,'ButtonDownFcn',@obj.LineClick);
            set(obj.LineHandle,'Color',obj.RULER_COLOR,'LineWidth',obj.LINE_WIDTH,'ButtonDownFcn',@obj.LineClick);
            
            %Create Labels at the side of the points
            obj.LabelHandle(1)=text(obj.COLORBAR_WIDTH-2,0.6,'');
            obj.LabelHandle(2)=text(obj.COLORBAR_WIDTH-2,obj.ColorMapSteps+0.5,'');
            
             %Setappdata which is used to prevent adding duplicate colorbars
            %or colorbars to colorbars
            setappdata(obj.ColorBarHandle,'ColorBarHandle',obj.ColorBarHandle);
            setappdata(obj.ColorBarHandle,'AxesHandle',obj.AxesHandle);
            setappdata(obj.ColorBarHandle,'ColorBarType','InteractiveColorbar');
            
            
            setappdata(obj.AxesHandle,'AxesHandle',obj.AxesHandle);
            setappdata(obj.AxesHandle,'ColorBarHandle',obj.ColorBarHandle);
            setappdata(obj.AxesHandle,'ColorBarType','AxesWithInteractiveColorbar');
            
            %Add context menu (right click to axes and image)
            obj.CreateContextMenu
            set(obj.ColorBarHandle,'UIContextMenu',obj.ContextMenu);
            set(Imageh,'UIContextMenu',obj.ContextMenu);
            set(obj.ColorBarHandle,'HandleVisibility','off');                      
        end
        
        function CreateContextMenu(obj)
            %Create context menu
            obj.ContextMenu=uicontextmenu('Callback',@SetCallBackMenu);
            for i=1:numel(obj.ColorMaps)
                uimenu(obj.ContextMenu,'Label',obj.ColorMaps{i},'CallBack',@obj.ContextMenuCallback);
            end
            
            ConnectHandle=uimenu(obj.ContextMenu,'Label','Connect Colorbar','CallBack',@obj.ContextMenuCallback);
            
            function SetCallBackMenu(~,~,~)
                if obj.UpdateConnectedColorBars;checkstatus='on';else checkstatus='off';end
                set(ConnectHandle,'Checked',checkstatus);
            end
        end
        
        %Update the point positions and colorbar image
        function Update(obj,UpdateConnections)
            if nargin==1;UpdateConnections=true;end
            %Calculate limits as fraction of Range
            Range=obj.CRange;
            Lim=obj.CLim;
            
            RangeSpan=Range(2)-Range(1);            
            LimFrac=(Lim-Range(1))/RangeSpan;
            
            LimFrac(2)=LimFrac(2);
            
            YSize=obj.ColorMapSteps;
            
            
            
            set(obj.PointHandle(1:2),'YData',YSize.*LimFrac(1));
            set(obj.PointHandle(3:4),'YData',YSize.*LimFrac(2));
            set(obj.LineHandle(1),'YData',[1 1].*YSize.*LimFrac(1));
            set(obj.LineHandle(2),'YData',[1 1].*YSize.*LimFrac(2));
            
            
            set(obj.ColorBarHandle,'CLim',LimFrac);
            set(obj.AxesHandle,'CLim',Lim);
            
            if (LimFrac(1)*obj.ColorMapSteps) >= (obj.ColorMapSteps-2)
                set(obj.PointHandle(3:4),'HitTest','off');
            else
                set(obj.PointHandle(3:4),'HitTest','on');
            end
            
            %Move Labels            
            set(obj.LabelHandle(1),'String',num2str(Lim(1),obj.NUMTEXT_FORMAT),'Position',[obj.COLORMAP_IMAGE_WIDTH+3, obj.ColorMapSteps.*LimFrac(1) 0]);
            set(obj.LabelHandle(2),'String',num2str(Lim(2),obj.NUMTEXT_FORMAT),'Position',[obj.COLORMAP_IMAGE_WIDTH+3, obj.ColorMapSteps.*LimFrac(2) 0]);
            
            
            if obj.UpdateConnectedColorBars && UpdateConnections                
                for i=1:numel(obj.ConnectedColorBarObjects)
                    if obj.ConnectedColorBarObjects(i).UpdateConnectedColorBars
                        set(obj.ConnectedColorBarObjects(i),'CRangeWithoutUpdateConnectedColorbarObjects',Range);
                        set(obj.ConnectedColorBarObjects(i),'CLimWithoutUpdateConnectedColorbarObjects',Lim);
                    end
                end
            end
            
            
        end
                                                              
        %Mouse Interactions
        function LineClick(obj,source,~)
            set(gcf,'WindowButtonMotionFcn',{@obj.MouseMove,source});
            set(gcf,'WindowButtonUpFcn',@obj.ButtonRelease);
        end
        
        function ButtonRelease(~,~,~)
            set(gcf,'WindowButtonMotionFcn','')
        end
        
        function MouseMove(obj,~,~,source)
            index=find(source==obj.LineHandle);
            if ~any(index);return;end
            
         
            newLim=obj.MouseCLim;
            
   
            
            switch index
                case {1}
                    obj.CLim(1)=newLim;
                    notify(obj,'UpdateByUser')
                case {2}
                    obj.CLim(2)=newLim;
                    notify(obj,'UpdateByUser')
            end
        end
                                        
        function ContextMenuCallback(obj,source,~)
            %Callback for context menu
            ColorMapStr=get(source,'Label');
            switch lower(ColorMapStr)                
                case 'colormap editor'
                    colormapeditor;
                case 'connect colorbar'
                    status=get(source,'Checked');
                    switch status
                        case 'on'
                            set(source,'Checked','off');
                            obj.UpdateConnectedColorBars=false;
                        case 'off'
                            set(source,'Checked','on');
                            obj.UpdateConnectedColorBars=true;
                            set(obj.ConnectedColorBarObjects,'CLim',obj.CLim);
                            %set(obj.ConnectedColorBarObjects,'CLim',obj.CLim);
                    end
                otherwise
                    obj.ColorMap=ColorMapStr;
            end
        end
                                                                                                                                                                          
    end
end