classdef transducer < handle
% Transducer class.
%   The transducer class is an implementation of various BK Medical 
%   transducers (8803, 8804, 8811, 8820e, 8823, 8862_1, 8862_2).
%    
%   Use: t = transducer(arg1, arg2,..., argN)
%
%   Input is specified in pairs of ('field name', 'value') i.e. 'id','8804'
%     Input can be 'test' for a test object.
%     
%   Output:
%     a transducer object:
%      * type: 'linear' / 'convex'
%      * name: i.e. '8804'
%      * nr_elements_x
%      * f0: center frequency           [Hz]
%      * bw: fractional bandwidth
%      * pitch                          [m]
%      * kerf                           [m]
%      * width                          [m]
%      * height                         [m]
%      * aclayerthickness               [m]
%      * elevation_focus                [m]
%      * element_positions  (x,y,z)     [m]
%      * element_directions (x,y,z)     [m]
%      * elevation_curvature_height     [m]       
%      * rcv_impulse_response_fs        [Hz]
%      * xmt_impulse_response_fs        [Hz]
%      * rcv_impulse_response - Default specified as the two-way IR if
%                               available, otherwise a delta function
%      * xmt_impulse_response - Default specified as a delta function
%      * ROC: Radius of Curvature (shell radius) [m]
%
%   Example of use:
%      t = transducer('test');
%
%   Methods:
%      * display: 
%         displays the object
%      * calculate_elevation_curvature_delay(c): 
%         calculates the time delay due to elevation focus
%         c is the speed of sound
%      * set_RCV_IR(IR,FS): 
%         set the rcv impulse response.
%      * set_XMT_IR(IR,FS): 
%         set the xmt impulse response.
%      * set_impulse_response_fs(fs): 
%         set the impulse response fs. Waveforms will be resampled 
%      * plot_IR:
%         plot the impulse response on current figure
%      * plot_elements:
%         plot the elements on current figure
%      * plot_GeoDebugger_elements:
%         plot the elements read from the GeoDebugger file on current figure
%      * ReadGeoDebugParameters:
%         Return the GeoDebugger parameters as a struct
%
% Martin Hemmsen (Ph.D)
% 18.03.11 (d.m.y)
% Ver. 03
% Changed the way the elements are placed in space. 
% 
    properties (SetAccess = 'public', GetAccess = 'public')
        id                      % Transducer ID        
        nr_elements_x           % Number of elements
        nr_elements_y           % Number of elements
        f0                      % Center frequency [Hz]
        bw                      % Fractional Bandwidth
        pitch                   % Pitch [m]
        kerf                    % Kerf [m]
        height                  % Height [m] (length in y-direction)
        aclayerthickness = 0    % Acoustic layer thickness [m]
    end % public properties

    properties (Dependent = true, SetAccess = 'private')
        % Width [m] (length in x-direction)   
        % Dependent on Pitch and Kerf
        % Width = Pitch - Kerf
        width    
        % Acoustic radius in [m]
        % AROC = ROC - aclayerthickness
        AROC
    end % Dependent private properties
    
    properties (SetAccess = 'public', GetAccess = 'public')
        elevation_focus         % Elevation focus [m]
        ROC                     % Outer radius of curvature (shell radius)
    end % public properties
 
    properties (SetAccess = 'private', GetAccess = 'public')
        rcv_impulse_response_fs % Receive IR sampling frequency
        xmt_impulse_response_fs % Transmit IR sampling frequency
        rcv_impulse_response    % Receive IR 
        xmt_impulse_response    % Transmit IR  
    end % private properties
   
    properties (Dependent = true, SetAccess = 'private')
       % Element positions (x,y,z)
       % Dependent on nr_elements_x, Pitch, AROC
       % First element at negative x-coordinates.
       % Last element at positive x-coordinates.
       % For more info see calc_element_position_LINEAR or
       % calc_element_position_CURVED
       element_positions       
       % Element directions 'Unit vector' (x,y,z)
       % Dependent on nr_elements_x, nr_elements_y, Pitch, AROC
       % Element direction is symmetric around 0.
       % Element direction increases from first to last element.
       % First element at negative x-coordinates.
       % Last element at positive x-coordinates.
       % For more info see calc_element_position_LINEAR or
       % calc_element_position_CURVED
       element_directions       
       % Element angles [rad]
       % Dependent on nr_elements_x, Pitch, AROC
       % Element angles is symmetric around 0.
       % Element angles increases from first to last element.
       % First element at negative x-coordinates.
       % Last element at positive x-coordinates.
       element_angles           
       % Element curvature height [m]
       elevation_curvature_height 
       % Transducer type 'linear' / 'convex'
       type
    end % Dependent private properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Methods                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods(Access='private',Hidden=true)
        function listfunctions(obj)
          [m n] = methods('transducer','-full');
%           [m n] = methods('transducer');
          for k=1:length(m)
            if (isempty(findstr(m{k},'L ')) && ...
                isempty(findstr(m{k},'TF ')) && ...
                isempty(findstr(m{k},'HM ')) && ...
                isempty(findstr(m{k},'Inherited')))
                fprintf('    ');
                fprintf(cell2mat(m(k)));
                fprintf('\n');
            end
          end
        end
    end
    
    methods(Access='private',Hidden=true)
        % FUNCTIONS DEFINED IN SEPARATE FILES
       [x, y, z, dir] = calc_element_position_LINEAR(obj,TotalLayerThickness,pitch,nr_elements_x,nr_elements_y);
       [x, z, dir] = calc_element_position_CURVED(obj,AROC,ROC,pitch,nr_elements_x);        
       parameter = ReadGeoDebugFile(obj,filename);
       [rect, center] = Create_Convex_Focused_Element(obj,ElementAngle,pitch,kerf,height,elevation_focus,ROC,ROCA,nr_sub_x,nr_sub_y);
       [rect] = Create_Linear_Focused_Element(obj,center,pitch,kerf,height,elevation_focus,aclayerthickness,nr_sub_x,nr_sub_y);        
    end
    
    methods
        % Construct an object
        function obj = transducer(varargin) 
            if nargin == 1 && strcmp(varargin{1}, 'test'), obj = transducer_test; return, end
            if nargin < 1, help(mfilename), error('Not enough input arguments.'), end

            % defaults     
            index = strcmpi(varargin,'id');
            index = find(index,1);
            if(~isempty(index))
                ID = varargin{index+1};
            else
                error('There must be an ID argument pair. i.e. ''id'',''8804'' ')
            end

            switch lower(ID) 
                case{'unspecified'}
                    obj.id                  = 'unspecified'; 
                    obj.nr_elements_x         = 0;                            % Number of elements in aperture 
                    obj.pitch               = 0;                            % Element pitch              [m]
                    obj.kerf                = 0;              	            % Element kerf               [m]
                    obj.height              = 0;                            % Element height             [m]
                    obj.elevation_focus     = 0;                            % Elevation focus (geometric)[m];       
                    obj.ROC                 = 0;                            % Outer radius of curvature    
                    obj.bw                  = 0;                            % Bandwidth [MHz]
                    obj.f0                  = 0;                            % Center frequency [MHz]

                    obj.rcv_impulse_response = 1;
                    obj.xmt_impulse_response = 1;
                    obj.rcv_impulse_response_fs = 40e6;
                    obj.xmt_impulse_response_fs = 40e6;
                    
                case{'8803'}
                    obj.id                  = '8803';
                    obj.nr_elements_x         = 128;                           % Number of elements in aperture 
                    obj.pitch               = 2.45e-004;                     % Element pitch              [m]
                    obj.kerf                = 0.035e-3;                  	 % Element kerf               [m]
                    obj.height              = 4e-3;                       	 % Element height             [m]
                    obj.elevation_focus     = 15 / 1000;                     % Elevation focus (geometric)[m];       
                    obj.ROC                 = 29.74 / 1000;                     % Outer radius of curvature    
                    obj.bw                  = 2.6e6;                         % Bandwidth [MHz]
                    obj.f0                  = 4e6;                           % Center frequency [MHz]
                    
                    obj.rcv_impulse_response = 1;
                    obj.xmt_impulse_response = 1;
                    obj.rcv_impulse_response_fs = 40e6;
                    obj.xmt_impulse_response_fs = 40e6;
                    
                case{'8804'}
                    obj.id                  = '8804';
                    obj.nr_elements_x         = 192;                           % Number of elements in aperture 
                    obj.pitch               = 2.08e-004;                     % Element pitch              [m]
                    obj.kerf                = 0.035  / 1000;              	 % Element kerf               [m]
                    obj.height              = 4.5 / 1000;               	 % Element height             [m]
                    obj.elevation_focus     = 20 / 1000;                     % Elevation focus (geometric)[m];       
                    obj.ROC                 = 0;                             % Outer radius of curvature    
                    obj.bw                  = 4.2e6;                         % Bandwidth [MHz]
                    obj.f0                  = 7e6;                           % Center frequency [MHz]
                                        
                    % Load 2-way impulse response
                    p = mfilename('fullpath');
                    i = strfind(p,filesep);
                    p = p(1:i(end));
                    load([p '8804' filesep 'impulse_bk8804.mat'],'pulse','fs') 
                    
                    % define rcv/xmt impulse response
                    obj.rcv_impulse_response = reshape(pulse(:),1,[]);
                    obj.xmt_impulse_response = 1;
                    obj.rcv_impulse_response_fs = fs;
                    obj.xmt_impulse_response_fs = fs;
                    
                case{'8811'}
                    obj.id                  = '8811';
                    obj.nr_elements_x         = 192;                           % Number of elements in aperture 
                    obj.pitch               = 0.26 / 1000;                   % Element pitch              [m]
                    obj.kerf                = 0.035 / 1000;              	 % Element kerf               [m]
                    obj.height              = 4  / 1000;                	 % Element height             [m]
                    obj.elevation_focus     = 15 / 1000;                     % Elevation focus (geometric)[m];       
                    obj.ROC                 = 0;                             % Outer radius of curvature    

                    obj.rcv_impulse_response = 1;
                    obj.xmt_impulse_response = 1;
                    obj.rcv_impulse_response_fs = 40e6;
                    obj.xmt_impulse_response_fs = 40e6;
                    
                case{'8670'}
                    obj.id                  = '8670';
                    obj.nr_elements_x         = 128;                           % Number of elements in aperture 
                    obj.pitch               = 0.3 / 1000;                    % Element pitch              [m]
                    obj.kerf                = 0.035 / 1000;              	 % Element kerf               [m]
                    obj.height              = 4  / 1000;                	 % Element height             [m]
                    obj.elevation_focus     = 20 / 1000;                     % Elevation focus (geometric)[m];       
                    obj.ROC                 = 0;                             % Outer radius of curvature    
                    obj.bw                  = 5.6e6;                         % Bandwidth [MHz]
                    obj.f0                  = 8e6;                           % Center frequency [MHz]
                                        
                    obj.rcv_impulse_response = 1;
                    obj.xmt_impulse_response = 1;
                    obj.rcv_impulse_response_fs = 40e6;
                    obj.xmt_impulse_response_fs = 40e6;
                    
                    % Load 2-way impulse response
                    p = mfilename('fullpath');
                    i = strfind(p,filesep);
                    p = p(1:i(end));
                    load([p '8670' filesep 'Unit_Step_Response.mat'],'pulse','fs') 
                    
                    % define rcv/xmt impulse response
                    obj.rcv_impulse_response = reshape(pulse(:),1,[]);
                    obj.xmt_impulse_response = 1;
                    obj.rcv_impulse_response_fs = fs;
                    obj.xmt_impulse_response_fs = fs;
                    
                case{'8820e'}
                    obj.id                  = '8820e';
                    obj.nr_elements_x         = 192;                           % Number of elements in aperture 
                    obj.pitch               = 0.33e-3;                       % Element pitch              [m]
                    obj.kerf                = 0.035e-3;              	     % Element kerf               [m]
                    obj.height              = 13 / 1000;               	     % Element height             [m]
                    obj.elevation_focus     = 65 / 1000;                     % Elevation focus (geometric)[m]  
                    obj.ROC                 = 61 / 1000;                     % (outer Radius Of Curvature)
                
                    % Load 2-way impulse response
                    p = mfilename('fullpath');
                    i = strfind(p,filesep);
                    p = p(1:i(end));
                    load([p '8820e' filesep 'Unit_Step_Response.mat'],'pulse','fs')
                    
                    % define rcv/xmt impulse response
                    obj.rcv_impulse_response = reshape(pulse(:),1,[]);
                    obj.xmt_impulse_response = 1;
                    obj.rcv_impulse_response_fs = fs;
                    obj.xmt_impulse_response_fs = fs;
                case{'8820e_sti'}
                    obj.id                  = '8820e_STI';
                    obj.nr_elements_x         = 192;                           % Number of elements in aperture 
                    obj.pitch               = 0.33e-3;                       % Element pitch              [m]
                    obj.kerf                = 0.035e-3;              	     % Element kerf               [m]
                    obj.height              = 13 / 1000;               	     % Element height             [m]
                    obj.elevation_focus     = 65 / 1000;                     % Elevation focus (geometric)[m]  
                    obj.ROC                 = 61377 / 1000;                     % (outer Radius Of Curvature)
                
                    % Load 2-way impulse response
                    p = mfilename('fullpath');
                    i = strfind(p,filesep);
                    p = p(1:i(end));
                    load([p '8820e' filesep 'Unit_Step_Response.mat'],'pulse','fs')
                    
                    % define rcv/xmt impulse response
                    obj.rcv_impulse_response = reshape(pulse(:),1,[]);
                    obj.xmt_impulse_response = 1;
                    obj.rcv_impulse_response_fs = fs;
                    obj.xmt_impulse_response_fs = fs;
                    
                case{'8823'}
                    obj.id                  = '8823';
                    obj.nr_elements_x         = 160;                           % Number of elements in aperture 
                    obj.pitch               = 0.196-3;                       % Element pitch              [m]
                    obj.kerf                = 0.035e-3;              	     % Element kerf               [m]
                    obj.height              = 12 / 1000;               	     % Element height             [m]
                    obj.elevation_focus     = 70 / 1000;                     % Elevation focus (geometric)[m]  
                    obj.ROC                 = 31 / 1000;                     % (outer Radius Of Curvature)
                
                    obj.rcv_impulse_response = 1;
                    obj.xmt_impulse_response = 1;
                    obj.rcv_impulse_response_fs = 40e6;
                    obj.xmt_impulse_response_fs = 40e6;
                    
                 case{'8862_1'}
                    obj.id                  = '8862_1';
                    obj.nr_elements_x         = 160;                           % Number of elements in aperture 
                    obj.pitch               = 0.16-3;                        % Element pitch              [m]
                    obj.kerf                = 0.050-3;              	     % Element kerf               [m]
                    obj.height              = 7 / 1000;               	     % Element height             [m]
                    obj.elevation_focus     = 40 / 1000;                     % Elevation focus (geometric)[m]  
                    obj.ROC                 = 25 / 1000;                     % (outer Radius Of Curvature)
                
                    obj.rcv_impulse_response = 1;
                    obj.xmt_impulse_response = 1;
                    obj.rcv_impulse_response_fs = 40e6;
                    obj.xmt_impulse_response_fs = 40e6;
                
                case{'8862_2'}
                    obj.id                  = '8862_2';
                    obj.nr_elements_x         = 160;                           % Number of elements in aperture 
                    obj.pitch               = 0.18-3;                        % Element pitch              [m]
                    obj.kerf                = 0.050-3;              	     % Element kerf               [m]
                    obj.height              = 7 / 1000;               	     % Element height             [m]
                    obj.elevation_focus     = 40 / 1000;                     % Elevation focus (geometric)[m]  
                    obj.ROC                 = 25 / 1000;                     % (outer Radius Of Curvature)
                
                    obj.rcv_impulse_response = 1;
                    obj.xmt_impulse_response = 1;
                    obj.rcv_impulse_response_fs = 40e6;
                    obj.xmt_impulse_response_fs = 40e6;
                        
                otherwise
                    fprintf('The transducer name was not recognized\n');
            end     
            % copy parameters from input to structure
            obj = va_arg(obj, varargin);
        end
        function GeoDebugger = ReadGeoDebugParameters(obj,varargin)
            switch(nargin)
                case 1
                    % Read GeoDebug parameter file
                    p = mfilename('fullpath');
                    i = strfind(p,filesep);
                    p = p(1:i(end));
                    if(exist([p obj.id filesep 'GeoDebug.txt']))
                        GeoDebugger = ReadGeoDebugFile(obj,[p obj.id filesep 'GeoDebug.txt']);
                    else
                       error('Parameters are not available - the GeoDebugger file was not found\n'); 
                    end
                case 2
                    % Read GeoDebug parameter file
                    p = mfilename('fullpath');
                    i = strfind(p,filesep);
                    p = p(1:i(end));
                    if(exist([p obj.id filesep varargin{1}]))
                        GeoDebugger = ReadGeoDebugFile(obj,[p obj.id filesep varargin{1}]);
                    else
                       error('Parameters are not available - the GeoDebugger file was not found\n'); 
                    end
                otherwise
                    error('wrong number of input arguments')
            end
            
        end
        function type = get.type(obj)
           % return type
           if(~isempty(obj.ROC))
               switch(obj.ROC)
                   case 0
                        type = 'linear';
                   otherwise
                       type = 'convex';
               end
           else
              error('type is unknown since ROC is not specified') 
           end
        end
        function width = get.width(obj)
           % calculate element width
           width = obj.pitch - obj.kerf;          % Element width              [m]
        end
        function AROC = get.AROC(obj)
            AROC = obj.ROC - obj.aclayerthickness;
        end
        function elevation_curvature_height = get.elevation_curvature_height(obj)
            % calculate elevation curvature height
            if(obj.elevation_focus == 0)
                elevation_curvature_height = 0;
            else
                elevation_curvature_height = obj.elevation_focus-cos(asin((obj.height/2)/obj.elevation_focus))*obj.elevation_focus;
            end
        end
        function value = calculate_elevation_curvature_delay(obj,c)
        % Calculate elevation curvature delay.
        % c: speed of sound
            value = obj.elevation_curvature_height/c;

        end        
        function element_directions = get.element_directions(obj)
           % calculate element directions
           switch(obj.type)
               case 'linear'
                   [dummyx, dummyy, dummyz, element_directions] = calc_element_position_LINEAR(obj,obj.aclayerthickness,obj.pitch,obj.nr_elements_x,obj.nr_elements_y);
               case 'convex'
                   [dummyx, dummyz, element_directions] = calc_element_position_CURVED(obj,obj.ROC-obj.aclayerthickness,obj.ROC,obj.pitch,obj.nr_elements_x);                
            end
            % First element at negative x-coordinates.
            % Last element at positive x-coordinates.
            % Element direction is symmetric around 0.
            % Element direction increases from first to last element.
            element_directions = element_directions(:);
            x = cos(element_directions);
            z = sin(element_directions);
            element_directions = [x(:),zeros(obj.nr_elements_x*obj.nr_elements_y,1),z(:)];
        end
        function element_angles = get.element_angles(obj)
           % calculate element directions
           switch(obj.type)
               case 'linear'
                   [dummyx, dummyy, dummyz, element_directions] = calc_element_position_LINEAR(obj,obj.aclayerthickness,obj.pitch,obj.nr_elements_x,obj.nr_elements_y);
               case 'convex'
                   [dummyx, dummyy, element_directions] = calc_element_position_CURVED(obj,obj.ROC-obj.aclayerthickness,obj.ROC,obj.pitch,obj.nr_elements_x);                
                   
            end
            % First element at negative x-coordinates.
            % Last element at positive x-coordinates.
            % Element direction is symmetric around 0.
            % Element direction increases from first to last element.
%             element_angles = element_directions(:) - pi/2;
            element_angles = element_directions(:);
        end      
        function element_positions = get.element_positions(obj)
           % calculate element positions
           switch(obj.type)
               case 'linear'
                   [x, y, z, dummydir] = calc_element_position_LINEAR(obj,obj.aclayerthickness,obj.pitch,obj.nr_elements_x,obj.nr_elements_y);
               case 'convex'
                   [x, z] = calc_element_position_CURVED(obj,obj.ROC-obj.aclayerthickness,obj.ROC,obj.pitch,obj.nr_elements_x);                
                   y = zeros(size(x));
           end
            
            % First nr_elements_x elements at positive y-coordinates.
            % First element at positive x-coordinates.
            % Last nr_elements_x elements at negative y-coordinates.
            % Last element at negative x-coordinates.
            % Element direction is symmetric around 0.
            % Element direction increases from first to last element.
            element_positions = [x(:),y(:),z(:)];
        end

        function plot_elements(obj,varargin)
            elements_to_plot = 1:obj.nr_elements_x;
            color = [0 0 1];
            displayname = 'Calculated element position';
            switch(nargin)
                case(2)
                    elements_to_plot = varargin{1};
                case(3)
                    elements_to_plot = varargin{1};
                    color = varargin{2};
                case(4)
                    elements_to_plot = varargin{1};
                    color = varargin{2};
                    displayname = varargin{3};
            end
                    
        % Displays the element positions on current figure
            pos = obj.element_positions;
            theta = obj.element_angles;
            % make the small line objects for each element
            a = cos(theta).*obj.width/2;
            b = sin(theta).*obj.width/2;

            x1 = pos(:,1)+b;
            x2 = pos(:,1)-b;

            y1 = pos(:,3)-a;
            y2 = pos(:,3)+a;
            
            legendAxListener = [];
            
            hold all
            for k = elements_to_plot
                if(k == elements_to_plot(1))
                    if(isempty(displayname))
                        try
                            legendListeners = get(gca,'ScribeLegendListeners');
                        catch
                            legend('-DynamicLegend');
                            legendListeners = get(gca,'ScribeLegendListeners');
                        end
                        legendAxListener = legendListeners.childadded;
                        set(legendAxListener,'Enable','off');
                    end
                    plot([x1(k) x2(k)]*1000,[y1(k) y2(k)]*1000,'color',color,'linewidth',6,...
                            'DisplayName',displayname);
                    try
                        legendListeners = get(gca,'ScribeLegendListeners');
                    catch
                        legend('-DynamicLegend');
                        legendListeners = get(gca,'ScribeLegendListeners');
                    end
                    legendAxListener = legendListeners.childadded;
                    set(legendAxListener,'Enable','off');
                else
                    plot([x1(k) x2(k)]*1000,[y1(k) y2(k)]*1000,'color',color,'linewidth',6);
                end
                    
            end
            % Re-enable the dynamic legend listener
            try
                set(legendAxListener,'Enable','on');
            end
            hold off
            
            axis ij
            if(isempty(get(get(gca,'xlabel'),'String')) && ...
               isempty(get(get(gca,'ylabel'),'String'))  && ...
               isempty(get(get(gca,'title'),'String')))     
               xlabel('Lateral position [mm]','fontsize',14)
               ylabel('Axial position [mm]','fontsize',14)
               title('Element positions','fontsize',16)
            end

        end       
        function plot_GeoDebugger_elements(obj)
        % Displays the GeoDebugger element positions on current figure
            theta = obj.element_angles;
            % make the small line objects for each element
            a = cos(theta).*obj.width/2;
            b = sin(theta).*obj.width/2;
            
            GeoDebugger = ReadGeoDebugParameters(obj);
            if(~isempty(GeoDebugger))
                pos = GeoDebugger.element_position_table(:,1:2);
                      
                x1 = pos(:,1)+b;
                x2 = pos(:,1)-b;

                y1 = pos(:,2)-a;
                y2 = pos(:,2)+a;

                hold all
                for k = 1:obj.nr_elements_x
                    if(k == 1)
                        plot([x1(k) x2(k)]*1000,[y1(k) y2(k)]*1000,'r','linewidth',3,...
                            'DisplayName','GeoDebugger element position');
                        try
                            legendListeners = get(gca,'ScribeLegendListeners');
                        catch
                            legend('-DynamicLegend');
                            legendListeners = get(gca,'ScribeLegendListeners');
                        end
                        legendAxListener = legendListeners.childadded;
                        set(legendAxListener,'Enable','off');
                    else
                        plot([x1(k) x2(k)]*1000,[y1(k) y2(k)]*1000,'r','linewidth',3);
                    end
                end
                hold off
                % Re-enable the dynamic legend listener
                try
                    set(legendAxListener,'Enable','on');
                end
            end
            
            axis ij
            if(isempty(get(get(gca,'xlabel'),'String')) && ...
               isempty(get(get(gca,'ylabel'),'String'))  && ...
               isempty(get(get(gca,'title'),'String')))     
               xlabel('Lateral position [mm]','fontsize',14)
               ylabel('Axial position [mm]','fontsize',14)
               title('Element positions','fontsize',16)
            end
        end       
        function plot_IR(obj)
        % Plot rcv/xmt impulse response. 
        % As default the rcv impulse response is the two-way IR and xmt IR
        % is a delta function.
            t_rcv = (0:size(obj.rcv_impulse_response,2)-1)*...
                1/obj.rcv_impulse_response_fs;
            
            NFFT = 1024*8;
            f_rcv = (0:NFFT-1)./NFFT*obj.rcv_impulse_response_fs-obj.rcv_impulse_response_fs/2;
            RCV_fs_IR = fftshift(abs(fft(obj.rcv_impulse_response,NFFT)));
            RCV_fs_IR = 20*log10(RCV_fs_IR./max(RCV_fs_IR(:)));
            
            figure
            subplot(2,1,1)
            hold on
            if(size(obj.rcv_impulse_response,2) > 1)
                plot(t_rcv*1e9,obj.rcv_impulse_response,'-bs','linewidth',2, 'MarkerEdgeColor','k',...
                        'MarkerFaceColor','r',...
                        'MarkerSize',4)
            else
                stem(t_rcv*1e9,obj.rcv_impulse_response,'linewidth',2)
            end

            subplot(2,1,2)
            hold on
            plot(f_rcv./1e6,RCV_fs_IR,'linewidth',2)
          
            
            t_xmt = (0:size(obj.xmt_impulse_response,2)-1)*...
                1/obj.xmt_impulse_response_fs;
            
            f_xmt = (0:NFFT-1)./NFFT*obj.xmt_impulse_response_fs-obj.xmt_impulse_response_fs/2;
            XMT_fs_IR = fftshift(abs(fft(obj.xmt_impulse_response,NFFT)));
            XMT_fs_IR = 20*log10(XMT_fs_IR./max(XMT_fs_IR(:)));
           
            subplot(2,1,1)
            if(size(obj.xmt_impulse_response,2) > 1)
                plot(t_xmt*1e9,obj.xmt_impulse_response,'-rs','linewidth',2, 'MarkerEdgeColor','k',...
                        'MarkerFaceColor','r',...
                        'MarkerSize',4)
            else
                stem(t_xmt*1e9,obj.xmt_impulse_response,'r','linewidth',2)
            end
            grid on
            xlabel('nsec','fontsize',14)
            ylabel('Amplitude','fontsize',14)
            title('Impulse response','fontsize',16)
            legend('rcv','xmt')
            axis tight
            hold off
            
            subplot(2,1,2)
            plot(f_xmt./1e6,XMT_fs_IR,'r','linewidth',2)
            grid on
            xlabel('MHz','fontsize',14)
            ylabel('amplitude dB','fontsize',14)
            title('Amplitude response','fontsize',16)
            legend('rcv','xmt')
            axis([0 max(obj.rcv_impulse_response_fs/2/1e6,obj.xmt_impulse_response_fs/2/1e6) -60 0])
            hold off
            
            drawnow
        end
        
        function set_RCV_IR(obj,IR,FS)
        % Set receive impulse response. 
        % IR: Impulse response - vector
        % FS: Sampling frequency
            if(nargin == 3)
                obj.rcv_impulse_response_fs = FS;
                obj.rcv_impulse_response = reshape(IR(:),1,[]);
            else
                error('Wrong number of arguments')
            end
        end
        function set_XMT_IR(obj,IR,FS)
        % Set transmit impulse response. 
        % IR: Impulse response - vector
        % FS: Sampling frequency
            if(nargin == 3)
                obj.xmt_impulse_response_fs = FS;
                obj.xmt_impulse_response = reshape(IR(:),1,[]);
            else
                error('Wrong number of arguments')
            end
        end      
        function set_impulse_response_fs(obj,val)
        % Set impulse response sample frequency  
        % rcv/xmt impulse response sampling frequency is set to the input
        % value. If the input value is different from the old sampling
        % frequency the rcv/xmt impulse response is resampled to match the
        % new sampling frequency.
            if(isempty(obj.xmt_impulse_response_fs))
                obj.xmt_impulse_response_fs = val;
            elseif(obj.xmt_impulse_response_fs ~= val)
                [interpolation_integer decimation_integer] = rat(val/obj.xmt_impulse_response_fs);
                if(length(obj.xmt_impulse_response) > 1)
                    % interpolate
                    pulse_fs = resample(obj.xmt_impulse_response,interpolation_integer,1);
                    % decimate
                    fir_order = min(floor(length(pulse_fs)/3)-1,13); % order > 13 is not recommended
                    pulse_fs = decimate(pulse_fs,decimation_integer,fir_order,'FIR');
                    obj.xmt_impulse_response = pulse_fs;
                end
                obj.xmt_impulse_response_fs = val;
%                 if(obj.display == 1)
%                     fprintf('Transducer class:\n')
%                     fprintf('  --> xmt_impulse_response_fs is set to %s.\n',num2str(obj.xmt_impulse_response_fs))
%                     fprintf('  --> XMT IR has been resampled.\n')
%                 end
            end
            
            if(isempty(obj.rcv_impulse_response_fs))
                obj.rcv_impulse_response_fs = val;
            elseif(obj.rcv_impulse_response_fs ~= val)
                if(length(obj.rcv_impulse_response) > 1)
                    % interpolate
                    pulse_fs = resample(obj.rcv_impulse_response,interpolation_integer,1);
                    % decimate
                    fir_order = min(floor(length(pulse_fs)/3)-1,13); % order > 13 is not recommended
                    pulse_fs = decimate(pulse_fs,decimation_integer,fir_order,'FIR');
                    obj.rcv_impulse_response = pulse_fs;
                end
                obj.rcv_impulse_response_fs = val;
%                 fprintf('Transducer class:\n')
%                 fprintf('  --> rcv_impulse_response_fs is set to %s.\n',num2str(obj.rcv_impulse_response_fs))
%                 fprintf('  --> RCV IR has been resampled.\n')                
            end
        end     
        
        function [rect,center] = create_aperture(obj,nr_sub_x,nr_sub_y)
            % Procedure for creating an aperture consisting of rectangles.
            % Output: rect: matrix of mathematical element description, see
            %               Field II manual
            %         center: matrix of center position of physical
            %               elements, see Field II manual
            rect = [];
            center = [];
            switch(obj.type)
                case{'linear'}
                    center = obj.element_positions;
                    for k = 1:obj.nr_elements_x
                        temp_rect = Create_Linear_Focused_Element(obj,center(k,:),...
                                                     obj.pitch,...
                                                     obj.kerf,...
                                                     obj.height,...
                                                     obj.elevation_focus,...
                                                     obj.aclayerthickness,...
                                                     nr_sub_x,...
                                                     nr_sub_y);
                        temp_rect(:,1) = k; % set element nr
                        rect = [rect; temp_rect]; % add element to matrix desciption
                    end
                case{'convex'}
                    for k = 1:obj.nr_elements_x
                        [temp_rect temp_center] = Create_Convex_Focused_Element(obj,obj.element_angles(k),...
                                                     obj.pitch,...
                                                     obj.kerf,...
                                                     obj.height,...
                                                     obj.elevation_focus,...
                                                     obj.ROC,...
                                                     obj.ROC-obj.aclayerthickness,...
                                                     nr_sub_x,...
                                                     nr_sub_y);
                        temp_rect(:,1) = k; % set element nr
                        rect = [rect; temp_rect]; % add element to matrix desciption
                        center = [center; temp_center];
                    end
                otherwise
                    error('wrong type')
            end
            

        end
        
        function show_aperture(obj,varargin)
            % Procedure for visualizing an aperture consisting of rectangles.
            % Input: If no input the aperture is visualized with 
            %        default nr_sub_x = 2 and nr_sub_y = 4.
            %        If input arg 1 is a matrix structure, from xdc_get
            %        (see Field II manual) this aperture is visualized
            %
            % Visualize rectangular elements
            switch(nargin)
                case 1 % no input arguments we just show the apperture
                       % with default nr_sub_x = 2 and nr_sub_y = 4
                    [rect,center] = create_aperture(obj,2,4);
                    apod_value = ones(size(rect,1),1);
                    corners = rect(:,2:13);
                case 2 % we asume input 1 is a structure from xdc_get
                    data = varargin{1};
                    apod_value = data(5,:)';
                    corners = data(11:22,:)';
            end
            
            colormap(cool(128));
            [N,M]=size(corners);
            % Do the actual display
            for i=1:N
                x=[corners(i,1), corners(i,10); corners(i,4), corners(i,7)]*1000;
                y=[corners(i,2), corners(i,11); corners(i,5), corners(i,8)]*1000;
                z=[corners(i,3), corners(i,12); corners(i,6), corners(i,9)]*1000;
                c=apod_value(i)*ones(2,2);
                hold on
%                s =  surf(x,y,z,c,'LineStyle','none');
               surf(x,y,z,c);
%                set(s,'EdgeColor',[0 0 1])
            end

            % Put som axis legends on
            Hc = colorbar;
            view(3)
            xlabel('x [mm]')
            ylabel('y [mm]')
            zlabel('z [mm]')
            grid
            axis('image')
            hold off
        end
        function display(obj) 
            % standard display class function
            fprintf('%s is an object of class %s:\n\n', inputname(1), class(obj))
            disp(obj)
            fprintf('  Methods:\n')
            obj.listfunctions();
            fprintf('\n');
        end
        function delete(obj)
            % destructor
        end  
    end % methods
end




%
% transducer_test()
%
function obj = transducer_test
  obj = transducer('id','8804');
end
