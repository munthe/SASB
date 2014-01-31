% ob = bft3_aperture(options) 
%
% Create an Aperture object by specifying a number of options as 
% string-arguments pairs.
%
% Options:
%   'type'         char,  'custom', 'linear_array',
%                         or 'convex_array' (default: 'custom')
%   'n_elements'   float, #elements
%   'pitch'        float, pitch
%   'kerf'         float, kerf   (only used for 'convex_array' 'type)
%   'width'        float, pitch  (only used for 'convex_array' 'type)
%   'radius'       float, radius (only used for 'convex_array' 'type)  
%   'pos'          float, positions of elements
%
% Properties:
%   pos is [#elements, 3]     float, set or get positions of the elements
%   focus is [1, 3]           float, set or get focus point (default: [])
%   center_focus is [1, 3]    float, set or get center focus (default:center of aperture)
%   delays is [1, #elements]  float, set or get receive delays (default: [])
%   ppwave is [1, 1]          bool, enable plane-wave calculation (default: false)
%
% Protected properties:
%   type             char,   type of aperture
%   Handle
%   Id
%
% Methods:
%   obj bft3_aperture(varargin)
%   out clone(obj)
%   delete(obj)
%   display(obj)
%   out focus_delays(obj),   get focus delays,   float[# elements, 1]
%
% If the option 'pos' is given, e.g.
%
% ah = bft3_aperture('pos',pos);
%
% 'n_elements', 'pitch' are ignored. To create a linear
% array facing, the xy-plane and positioned centrally and aligned with
% the x-axis, set the option 'type' to 'linear_array' together with the
% options 'pitch', and 'n_elements'. To create a convex array
% facing the xy-plane and positioned centrally and aligned with the
% x-axis, set the option 'type' to 'convex_array' together with the
% options 'width', 'kerf', 'radius', and 'n_elements'.
%
% Example:
%   fs = 30e6; f0 = 5e6; c = 1480; lambda = c/f0;
%   pitch = lambda;
%   globals = bft3_system('fs',fs,'c',c);
%   ah = bft3_aperture('type','linear_array','pitch',pitch,...
%                      'n_elements',64);
%
% $Id: bft3_aperture.m,v 1.64 2011-08-30 20:05:31 jmh Exp $
%

% @file bft3_aperture.m
% @brief Aperture class
%>
%>
%======================================================================
%> @brief Aperture class
%>
%> Create an Aperture object by specifying a number of options as 
%> string-arguments pairs.
%>
%> @par Options:\n
%> <table rules="none">
%>  <tr><td>@b 'type'</td><td>char</td><td>'custom' or 'linear_array' (default: 'custom')</td></tr>
%>  <tr><td>@b 'n_elements'</td><td>float</td><td># elements</td></tr>
%>  <tr><td>@b 'pitch'</td><td>float</td><td>pitch</td></tr>
%>  <tr><td>@b 'width'</td><td>float</td><td>width (only specified
%>for 'convex_array')</td></tr>
%>  <tr><td>@b 'kerf'</td><td>float</td><td>kerf (only specified
%>for 'convex_array')</td></tr>
%>  <tr><td>@b 'radius'</td><td>float</td><td>radius (only specified
%>for 'convex_array')</td></tr>
%>  <tr><td>@b 'pos'</td><td>float</td><td>positions of elements</td></tr>
%> </table>
%>
%> @par Properties:\n
%> <table rules="none">
%>  <tr><td>@ref pos </td><td>float[# elements, 3]</td><td>set or get positions of the elements                 </td></tr>
%>  <tr><td>@ref focus </td><td>float[1, 3]</td><td>set or get focus point (default: [])                 </td></tr>
%>  <tr><td>@ref center_focus</td><td>float[1, 3]</td><td>set or get center focus (default:center of aperture) </td></tr>
%>  <tr><td>@ref delays</td><td>float[1, # elements]</td><td>set or get receive delays (default: [])              </td></tr>
%>  <tr><td>@ref ppwave</td><td>bool</td><td>enable or disable plane-wave calculation (default: false)</td></tr>
%> </table>
%>
%> @par Read-only properties:\n
%> <table rules="none">
%>  <tr><td>@ref type   </td><td>char </td><td>type of aperture</td></tr>
%>  <tr><td>@ref Handle </td><td>uint(32/64)</td><td>pointer</td></tr>
%>  <tr><td>@ref Id     </td><td>uint(32/64)</td><td>unique identifier</td></tr>
%> </table>
%>
%> @par Methods:\n
%>   obj @ref bft3_aperture (varargin)\n
%>   out @ref clone (obj)\n
%>   @ref delete (obj)\n
%>   @ref display (obj)\n
%>   out @ref focus_delays (obj),   get focus delays,   float[# elements, 1]
%>
%>@par Example:\n
%>If the option 'pos' is given, e.g.
%>@code
%>ah = bft3_aperture('pos',pos);
%>@endcode
%>'n_elements', and 'pitch' are ignored.\n\n To create a linear
%>array facing, the xy-plane and positioned centrally and aligned with
%>the x-axis, set the option 'type' to 'linear_array' together with the
%>options 'pitch', and 'n_elements'.
%>@code
%>fs = 30e6; f0 = 5e6; c = 1480; lambda = c/f0;
%>pitch = lambda;
%>globals = bft3_system('fs',fs,'c',c);
%>ah = bft3_aperture('type','linear_array','pitch',pitch,...
%>                   'n_elements',64);
%>@endcode
%>To create a convex array
%>facing the xy-plane and positioned centrally and aligned with the
%>x-axis, set the option 'type' to 'convex_array' together with the
%>options 'width', 'kerf', 'radius', and 'n_elements'.
%>@code
%>fs = 30e6; f0 = 5e6; c = 1480; lambda = c/f0;
%>pitch = lambda; kerf = pitch/10; width = pitch-kerf;
%>radius = 20/1000;
%>globals = bft3_system('fs',fs,'c',c);
%>ah = bft3_aperture('type','convex_array','width',width,...
%>                   'kerf',kerf,'radius',radius,'n_elements',64);
%>@endcode
%>  $Id: bft3_aperture.m,v 1.64 2011-08-30 20:05:31 jmh Exp $
%>
% ======================================================================
classdef bft3_aperture < handle
  %
  % TODO:
  %  - bft3_va_arg should notify if properties have been changed, e.g.
  %    fs is the same as default
  %  - type should be contained in C++
  %

  % Constant properties
  properties (Constant, GetAccess = 'private', SetAccess = ...
              'private')
    %> @brief Library execution unit
    %>
    %> Name of mex-file called by the @ref bft3_aperture class
      
    % Library execution unit
    mexname = 'bft3_mex';%'aperture_mex';
  end % properties
    
  % Read-only properties    
  properties (SetAccess = 'private', GetAccess = 'public')
    %> @brief uint(32/64) object handle (read-only)

    % uint(32/64) object handle (read-only)
    Handle
    %> @brief uint(32/64) unique identifier of data content
    %>(read-only)
    
    % uint(32/64) unique identifier of data content (read-only)
    Id
  end % properties

  properties (Dependent = true, SetAccess = 'private', GetAccess = 'public')
    %> @brief char string classifier (read-only)
    %>
    %> String classifier set if aperture is constructed by
    %> specifying a 'type' together with a number of options.
    %> Possible apertures types include 'linear_array',
    %> 'convex_array', and 'custom'.
      
    % char string classifier (read-only)
    %
    % String classifier set if aperture is constructed by
    % specifying a 'type' together with a number of options
    % Possible apertures types include 'linear_array',
    % 'convex_array', and 'custom'.
    type          
  end % properties
  
  % Dependent properties
  properties (Dependent = true, SetAccess = 'public', GetAccess = 'public')
    %> @brief float[# elements 3] positions

    % float[# elements 3] positions
    pos           
    %> @brief float[1 3] focus or virtual source
    %>
    %> Set this property to a point in space float[1 3] for introducing a virtual source for
    %> the aperture. Virtual source can be introduced for the
    %> transmit or the receive aperture or both. When no virtual
    %> source is present the value should be empty
      
    % float[1 3] focus or virtual source
    % Set this property to a point in space float[1 3] for introducing a virtual source for
    % the aperture. Virtual source can be introduced for the
    % transmit or the receive aperture or both. When no virtual
    % source is present the value should be empty
    %
    focus
    %> @brief float[1 3] reference position for time-of-flight (TOF)
    %> calculations.
    %>
    %> Set this to the position corresponding to the sample at time
    %> equal to zero. For Field II simulations, this corresponds to
    %> the position you have set using xdc_set_center_focus.

    % float[1 3] reference position for time-of-flight (TOF) calculations.
    %
    % Set this to the position corresponding to the sample at time
    % equal to zero. For Field II simulations, this corresponds to
    % the position you have set using xdc_set_center_focus.
    center_focus  
    %> @brief float[1 3] orientation used for transmit apodization
    %> or plane-waves

    % float[1 3] orientation used for transmit apodization
    % or plane-waves 
    orientation
    %> @brief float[1 # elements] time delays
    %>
    %> Constant values added to any delays

    % float[1 # elements] time delays
    %
    % Constant values added to any delays
    delays
    %> @brief bool Plane-wave option
    %>
    %> Enable this property for plane-wave
    %> time-of-flight (TOF) calculation. Time-of-flight is computed
    %> using the distance from the plane containing @ref
    %> center_focus and perpendicular to
    %> @ref bft3_aperture::orientation of the aperture to
    %> the focal points. If a virtual source is used, i.e.
    %> @ref focus is non-zero, the time-of-flight is computed like
    %> in the case of a virtual source, but the distance from the
    %> virtual source to the focal points is now instead computed
    %> using the distance from the plane containing @ref focus
    %> to the the focal points.
      
    % Plane-wave option
    %
    % Enable this property for plane-wave
    % time-of-flight (TOF) calculation. Time-of-flight is computed
    % using the distance from the plane containing
    % center_focus and perpendicular to
    % bft3_aperture::orientation of the aperture to
    % the focal points. If a virtual source is used, i.e.
    % focus is non-zero, the time-of-flight is computed like
    % in the case of a virtual source, but the distance from the
    % virtual source to the focal points is now instead computed
    % using the distance from the plane containing focus
    % to the the focal points.
    ppwave
  end % properties

  properties (Dependent = true, SetAccess = 'private', GetAccess = 'private')
    %> Sampling frequency
      
    % Sampling frequency
    fs
    %> Speed of sound
      
    % Speed of sound
    c
    %> Center frequency

    % Center frequency
    f0            
  end % properties
  
  methods
    % ======================================================================
    %> @brief Class constructor
    %>
    %> Create an Aperture object by specifying a number of options as 
    %> string-argument pairs.
    %>
    %> @par Options:
    %> <table rules="none">
    %>  <tr><td>@b 'type'</td><td>char</td><td>'custom' or 'linear_array' (default: 'custom')</td></tr>
    %>  <tr><td>@b 'n_elements'</td><td>float</td><td># elements</td></tr>
    %>  <tr><td>@b 'pitch'</td><td>float</td><td>pitch</td></tr>
    %>  <tr><td>@b 'width'</td><td>float</td><td>width (only specified
    %>for 'convex_array')</td></tr>
    %>  <tr><td>@b 'kerf'</td><td>float</td><td>kerf (only specified
    %>for 'convex_array')</td></tr>
    %>  <tr><td>@b 'radius'</td><td>float</td><td>radius (only specified
    %>for 'convex_array')</td></tr>
    %>  <tr><td>@b 'pos'</td><td>float</td><td>positions of elements</td></tr>
    %> </table>
    %>
    %> @param varargin a number of options as string-argument pairs
    %>
    %> @return instance of the bft3_aperture class.
    %>@par Example:
    %>If the option 'pos' is given, e.g.
    %>@code
    %>ah = bft3_aperture('pos',pos);
    %>@endcode
    %>'n_elements', and 'pitch' are ignored.\n\n To create a linear
    %>array facing, the xy-plane and positioned centrally and aligned with
    %>the x-axis, set the option 'type' to 'linear_array' together with the
    %>options 'pitch', and 'n_elements'.
    %>@code
    %>fs = 30e6; f0 = 5e6; c = 1480; lambda = c/f0;
    %>pitch = lambda;
    %>globals = bft3_system('fs',fs,'c',c);
    %>ah = bft3_aperture('type','linear_array','pitch',pitch,...
    %>                   'n_elements',64);
    %>@endcode
    %>To create a convex array
    %>facing the xy-plane and positioned centrally and aligned with the
    %>x-axis, set the option 'type' to 'convex_array' together with the
    %>options 'width', 'kerf', 'radius', and 'n_elements'.
    %>@code
    %>fs = 30e6; f0 = 5e6; c = 1480; lambda = c/f0;
    %>pitch = lambda; kerf = pitch/10; width = pitch-kerf;
    %>radius = 20/1000;
    %>globals = bft3_system('fs',fs,'c',c);
    %>ah = bft3_aperture('type','convex_array','width',width,...
    %>                   'kerf',kerf,'radius',radius,'n_elements',64);
    %>@endcode
    % ======================================================================
    function obj = bft3_aperture(varargin)
      % Create an bft3_aperture object by specifying a number of options as 
      % string-argument pairs.
      %
      %  Calling:  ah = bft3_aperture(options);
      %
      %  Options:
      %    'type'           char,  'linear_array', 'convex_array', or (default:'custom') 
      %    'n_elements'     float, # elements
      %    'pitch'          float[m], pitch
      %    'width'           float, (only specified for 'convex_array')        
      %    'kerf'           float, (only specified for 'convex_array')
      %    'radius'         float, (only specified for 'convex_array')
      %    'pos'            float[m], positions of elements
        
      % Defaults
      st.pos         = [1 2 3 ; 4 5 6 ; 7 8 9 ; 10 11 12];
      st.type        = 'custom';
      st.n_elements  = 64;
      st.pitch       = 0.001;
      st.height      = 0.005;
      st.radius      = 20/1000;
      st.kerf        = 0.05/1000;
      st.width       = 0.220/1000;
      st.c           = 1540;
      st.fs          = 70e6;
      st.ppwave      = false;
      st.delays      = [];
      st.focus       = [];
      st.name        = [];
      st.f0          = 3e6;
      
      if (nargin == 1)
        if (strcmp(varargin{1},'test'))
          % Unit test of aperture
          obj = bft3_aperture_test();
          return;
        elseif isstruct(varargin{1})
          aptype = varargin{1}.type;
          st = varargin{1};
          st.type = 'custom';
        elseif (all(size(varargin{1})==[1 1]) &&...
                (strcmp(class(varargin{1}),'uint32') || ...
                  strcmp(class(varargin{1}),'uint64')))
          % Aperture returned by apodization object
          obj.Handle = varargin{1};
          return;
        else
          help(mfilename), error(mfilename)
        end
      else
        st = bft3_va_arg(st,varargin);
      end
      
      if ~isempty(st.name)
        if strcmp(st.name,'8820e')
          st.type = 'convex_array';
          st.radius     = 0.0603;
          st.kerf       = 3.5e-5;
          st.pitch      = 3.3e-4;
          st.n_elements = 192;
          st.f0         = 3.5e6;
        elseif strcmp(st.name,'stv4')
          st.type       = 'linear_array';
          st.n_elements = 128;
          st.pitch      = 2.2e-4;
          st.f0         = 3.5e6;
        elseif strcmp(st.name,'8804')
          st.type = 'linear_array';
          st.n_elements = 192;
          st.pitch  = 0.208/1000;
          st.f0     = 7e6;          
        elseif strcmp(st.name,'stv7')
          st.type = 'linear_array';          
          st.n_elements = 128;
          st.pitch  = 1.1e-4;
          st.f0     = 7e6;
        end
      end
      
      if strcmp(st.type,'custom')
        % Manual aperture
        pos = st.pos;
        eval(['obj.Handle=', bft3_aperture.mexname,...
              '(''aperture,ctor,manual'',bft3_conv_float(pos));']);
      elseif  strcmp(st.type,'linear_array')
        % Field II aperture
        eval(['obj.Handle=', bft3_aperture.mexname,...
              '(''aperture,ctor,field'',',...
              'st.type, int32(st.n_elements),'...
              'bft3_conv_float(st.pitch),' ...
              'bft3_conv_float(st.f0));']);
      elseif  strcmp(st.type,'convex_array')
        % Field II aperture
        eval(['obj.Handle=', bft3_aperture.mexname,...
              '(''aperture,ctor,field,convex'',',...
              'st.type, int32(st.n_elements),'...
              'bft3_conv_float(st.width),' ...
              'bft3_conv_float(st.kerf),' ...
              'bft3_conv_float(st.radius),' ...
              'bft3_conv_float(st.f0));']);
      else
        error('Unsupported aperture type')
      end
      
      % Update when struct is used for creation
      if isstruct(varargin{1})
        obj.focus        = st.focus;
        obj.center_focus = st.center_focus;
        obj.delays       = st.delays;
        obj.fs           = st.fs;
        obj.c            = st.c;
        obj.type         = aptype;
        obj.ppwave       = st.ppwave;
      else
        obj.focus        = st.focus;
      %        obj.center_focus = st.center_focus;
        obj.delays       = st.delays;
        obj.ppwave       = st.ppwave;
      end
    end
  
    % ======================================================================
    %> @brief Display function
    %>
    %> Display the properties and member functions of the class
    %>    
    %> @param obj instance of the bft3_aperture class.
    % ======================================================================
    function display(obj)
      % Display
      name = inputname(1);
      fprintf('%s is an object of class %s:\n\n', name, ...
              class(obj))
      disp(obj)
      fprintf('  Methods:\n')
      bft3_listfunctions(obj);
      fprintf('\n');
    end
  end % methods
  methods (Access = 'private')
    % ======================================================================
    %> @brief Clone aperture
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @retval (deep) copy
    % ======================================================================
    function out = clone(obj)
      % Clone aperture
      eval(['Handle=', bft3_aperture.mexname,...
            '(''aperture,get,clone'',obj.Handle);']);
      out = bft3_aperture(Handle);
    end
  end % methods
  methods
    % ======================================================================
    %> @brief Element positions
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @retval out float[\#elements 3]
    % ======================================================================
    function out = get.pos(obj)
      % Element positions
      eval(['out=', bft3_aperture.mexname,...
            '(''aperture,get,pos'',obj.Handle);']);
    end
  
    % ======================================================================
    %> @brief Speed of sound
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @retval out float speed of sound
    % ======================================================================
    function out = get.c(obj)
      % Speed of sound
      eval(['out=', bft3_aperture.mexname,...
            '(''aperture,get,c'',obj.Handle);']);
    end

    % ======================================================================
    %> @brief Center focus
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @retval out center focus
    % ======================================================================
    function out = get.center_focus(obj)
      % Center focus
      eval(['out=', bft3_aperture.mexname,...
            '(''aperture,get,center_focus'',obj.Handle);']);
    end
    
    % ======================================================================
    %> @brief Orientation
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @retval out orientation
    % ======================================================================
    function out = get.orientation(obj)
      % Center focus
      eval(['out=', bft3_aperture.mexname,...
            '(''aperture,get,orientation'',obj.Handle);']);
    end

    % ======================================================================
    %> @brief Delays
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @retval out delays
    % ======================================================================
    function out = get.delays(obj)
      % Delays
      eval(['out=', bft3_aperture.mexname,...
            '(''aperture,get,delays'',obj.Handle);']);
    end

    % ======================================================================
    %> @brief Focus
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @retval out float[1 3] position of focus or virtual source
    % ======================================================================
    function out = get.focus(obj)
      % Focus
      eval(['out=', bft3_aperture.mexname,...
            '(''aperture,get,focus'',obj.Handle);']);
    end
  
    % ======================================================================
    %> @brief Center frequency
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @retval out float center frequency
    % ======================================================================
    function out = get.f0(obj)
      % Center frequency
      eval(['out=', bft3_aperture.mexname,...
            '(''aperture,get,f0'',obj.Handle);']);
    end
  
    % ======================================================================
    %> @brief Sampling frequency
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @retval out float sampling frequency
    % ======================================================================
    function out = get.fs(obj)
      % Sampling frequency
      eval(['out=', bft3_aperture.mexname,...
            '(''aperture,get,fs'',obj.Handle);']);
    end

    % ======================================================================
    %> @brief Id
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @retval out uint(32/64) unique identifier of aperture contents
    % ======================================================================
    function out = get.Id(obj)
      % Id
      eval(['out=', bft3_aperture.mexname,...
            '(''aperture,get,id'',obj.Handle);']);
    end

    % ======================================================================
    %> @brief Type
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @retval out char type of aperture
    % ======================================================================
    function out = get.type(obj)
      % Type
      eval(['out=', bft3_aperture.mexname,...
            '(''aperture,get,type'',obj.Handle);']);
    end

    % ======================================================================
    %> @brief Plane-parallel wave
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @retval out bool
    % ======================================================================
    function out = get.ppwave(obj)
      % Plane-parallel wave
      eval(['out=', bft3_aperture.mexname,...
            '(''aperture,get,ppwave'',obj.Handle);']);
    end
    
    % ======================================================================
    %> @brief Set positions
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @param data float[\#elements 3]
    % ======================================================================
    function set.pos(obj,data)
      % Positions
      eval([bft3_aperture.mexname,...
            '(''aperture,set,pos'',obj.Handle,bft3_conv_float(data));']);
    end
  
    % ======================================================================
    %> @brief Set center focus
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @param data float[1 3]
    % ======================================================================
    function set.center_focus(obj,data)
      % Center focus
      eval([bft3_aperture.mexname,...
            '(''aperture,set,center_focus''',...
            ',obj.Handle,bft3_conv_float(data));']);
    end
    
    % ======================================================================
    %> @brief Set orientation
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @param data float[1 3]
    % ======================================================================
    function set.orientation(obj,data)
      % Orientation
      eval([bft3_aperture.mexname,...
            '(''aperture,set,orientation''',...
            ',obj.Handle,bft3_conv_float(data));']);
    end

    % ======================================================================
    %> @brief Set delays
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @param data [1 \#elements]
    % ======================================================================
    function set.delays(obj,data)
      % Delays
      eval([bft3_aperture.mexname,...
            '(''aperture,set,delays''',...
            ',obj.Handle,bft3_conv_float(data));']);
    end
    
    % ======================================================================
    %> @brief Set focus
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @param data float[1 3]
    % ======================================================================
    function set.focus(obj,data)
      % Focus
      eval([bft3_aperture.mexname,...
            '(''aperture,set,focus''',...
            ',obj.Handle,bft3_conv_float(data));']);
    end
  end % methods
  methods
    % ======================================================================
    %> @brief Set sampling frequency
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @param data float
    % ======================================================================
    function set.fs(obj,data)
      % Set sampling frequency
      eval([bft3_aperture.mexname,...
            '(''aperture,set,fs''',...
            ',obj.Handle,bft3_conv_float(data));']);
    end
    
    % ======================================================================
    %> @brief Set speed of sound
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @param data float
    % ======================================================================
    %
    function set.c(obj,data)
      % Set speed of sound
      eval([bft3_aperture.mexname,...
            '(''aperture,set,c''',...
            ',obj.Handle,bft3_conv_float(data));']);
    end

    % ======================================================================
    %> @brief Set center frequency
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @param data float
    % ======================================================================
    function set.f0(obj,data)
      % Set center frequency
      eval([bft3_aperture.mexname,...
            '(''aperture,set,f0''',...
            ',obj.Handle,bft3_conv_float(data));']);
    end

    % ======================================================================
    %> @brief Set type (read-only)
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @param data char
    % ======================================================================
    function set.type(obj,data)
      eval([bft3_aperture.mexname,...
            '(''aperture,set,type''',...
            ',obj.Handle, data);']);
    end
    
  end
  methods
    % ======================================================================
    %> @brief Set ppwave option
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @param data bool
    % ======================================================================
    function set.ppwave(obj,data)
      eval([bft3_aperture.mexname,...
            '(''aperture,set,ppwave'',obj.Handle,data);']);
    end
    
    % ======================================================================
    %> @brief Retrieve delay values for the current focus set using
    %> the @ref focus property
    %>
    %> @param obj instance of the bft3_aperture class.
    %> @retval focus delays corresponding to the current focus
    % ======================================================================
    function out = focus_delays(obj)
      % Get focus delays
      eval(['out=', bft3_aperture.mexname,...
            '(''aperture,get,focus_delays'',obj.Handle);']);
    end

    % ======================================================================
    %> @brief Class destructor
    %>
    %> Delete method are called before an object of the class is destroyed 
    %> @param obj instance of the bft3_aperture class.
    % ======================================================================
    function delete(obj)
      % Class destructor
      try
        eval([bft3_aperture.mexname,...
            '(''aperture,dtor'',obj.Handle);']);
      catch me
        % Throws warning if garbage collector for dll has cleaned the object
        bft3_warn(['Garbage collector has already freed the memory:\n',...
                   me.message])
      end
    end
  end  % methods
end % class

% ======================================================================
%> @brief Unit test of the bft_aperture class
%>
%> Function included for testing consistency. Throws an error in
%> case of in-consistency
%>
%> @return obj instance of the bft3_aperture class.
% ======================================================================
function obj = bft3_aperture_test()
  % Unit test of the bft_aperture class
  c = 1540;
  f0 = 5e6;
  lambda = c/f0;
  pitch = lambda;
  n_elements = 64;
  obj = bft3_aperture('type','linear_array','pitch',pitch,...
                      'n_elements',n_elements);

  wx = (n_elements-1)/2;
  refpos = ((0:n_elements-1) - wx)*pitch;
  if (sum(abs(refpos-obj.pos(:,1)')) > eps)
    bft3_warn('linear_array element positions are wrong')
  end
  
  % Change number of emissions by specifying orientation of larger
  % dimension. Verify that the already defined center foci are
  % unchanged
  obj.center_focus = [1 1 1];
  obj.orientation  = 2*ones(2,3);
  if (~all(size(obj.center_focus) == [2 3]) | ...
      ~all(obj.center_focus == [1 1 1 ; 0 0 0]))
    error('Unit test failed');
  end
  obj.center_focus = [0 0 0];
  obj.orientation  = [0 0 0];
end
%> <table rules="none">
%>  <tr><td> </td><td> </td><td> </td></tr>
%>  <tr><td> </td><td> </td><td> </td></tr>
%>  <tr><td> </td><td> </td><td> </td></tr>
%>  <tr><td> </td><td> </td><td> </td></tr>
%> </table>
