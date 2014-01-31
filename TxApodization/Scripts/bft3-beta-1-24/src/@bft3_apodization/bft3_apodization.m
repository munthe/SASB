% ob = bf3_apodization(aperture, ref, distances, values)
%
% Create an apodization object using four arguments or by specifying
% a number of options as string-arguments pairs.
%
% Calling: ob = bf3_apodization(aperture, ref, distances, values)        
%
% Parameters:
%    aperture                        class bft3_aperture
%    ref is [1 3]                    float, apodization reference point
%    distances is [1, # dist]        float, distances to reference point
%    values is [# elements, # dist]  float, apodization values
%
%  or
%          ob = bft3_apodization(aperture, options)
% 
% Options:
%    ref is [1,3]                   float, Apodization reference point
%    distances is [1, # dist]       float, Distances to this reference
%    values is [# elements, # dist] float, Apodization values (window functions)
%
% Properties: (can be modified after construction)
%    ref is [1,3]                 float, apodization reference point
%    parametric                   bool, enable parametric apodization,
%                                 i.e enable apodization windows
%                                 each specified by a distance and
%                                 a window. Each windows are active when
%                                 we are further away than the
%                                 distance.
%    distances is [1, #dist]      float, distances to reference point
%    values is [#elements, #dist] float, apodization values (window
%                                 functions)
%    fixed                        bool, enable fixed apodization.
%                                 Apodization window of width
%                                 n_active_elements and type
%                                 'window' with an orientation
%                                 specified by 3 Euler angles (only
%                                 supported for bft3_sampled_image)
%    n_active_elements            float, number of active elements
%    dynamic                      bool, enable dynamic apodization,
%                                 specified by an F-number and a
%                                 window type.
%    f                            float, F-number used for dynamic apodization
%    window                       char, window function used for
%                                 dynamic and fixed apodization,
%                                 'Rectwin', 'Hamming', 'Hann',
%                                 'Blackman', 'Tukey', 'Gaussian',
%                                 or 'Bartlett'.
%    window_parameter             Some windows require a
%                                 parameter. For the 'Gaussian'
%                                 window, this is the inverse
%                                 std. deviation. For the 'Tukey'
%                                 window, this is ratio of taper to
%                                 constant sections normalized to
%                                 (0,1); 0 (Hanning), 1 (Rectwin)
%
% Read-only properties:
%    Handle
%    aperture
%    Id
%
% Methods:
%    obj bft3_apodization(varargin)
%    obj clone(obj)
%    delete(obj)
%    display(obj)
%
% $Id: bft3_apodization.m,v 1.59 2012-01-19 10:59:54 jmh Exp $
%
  
% @file bft3_apodization.m
% @brief Apodization class 
% ======================================================================
%> @brief Apodization class
%>
%> Create an apodization object using four arguments or by specifying
%> a number of options as string-arguments pairs.
%>
%> @par Calling:
%> \n ob = bf3_apodization(aperture, ref, distances, values)        
%>
%>  or
%> @par
%>          ob = bft3_apodization(aperture, options)
%>
%> @par Parameters:
%> <table rules="none">
%>  <tr><td>aperture                         </td><td>class bft3_aperture</td><td></td></tr>
%>  <tr><td>ref is [1 3]                     </td><td>float</td><td>apodization reference point</td></tr>                             
%>  <tr><td>distances is [1, \#dist]         </td><td>float</td><td>distances to reference point</td></tr>
%>  <tr><td>values is [\#elements, \#dist]   </td><td>float</td><td>apodization values</td></tr>
%> </table>                                                         
%>      
%> @par Options:
%> <table rules="none">
%>  <tr><td>ref is [1,3]                    </td><td>float </td><td>Apodization reference point           </td></tr>
%>  <tr><td>distances is [1, \#dist]        </td><td>float </td><td>Distances to this reference           </td></tr>
%>  <tr><td>values is [\#elements, \#dist]  </td><td>float </td><td>Apodization values (window functions) </td></tr>
%> </table>
%>
%> @par Properties: (can be modified after construction)
%> <table rules="none">
%>  <tr><td>ref is [1,3]</td><td>float</td><td>apodization reference point</td></tr>
%>  <tr><td>parametric</td><td>bool</td><td>enable parametric apodization,
%>                                 i.e enable apodization windows
%>                                 each specified by a distance and
%>                                 a window. Each windows are active when
%>                                 we are further away than the
%>                                 distance.</td></tr>
%>  <tr><td>distances is [1, \#dist]</td><td>float</td><td>distances to reference point</td></tr>
%>  <tr><td>values is [\#elements, \#dist]</td><td>float</td><td>apodization values (window
%>                                 functions)</td></tr>
%>  <tr><td>fixed</td><td>bool</td><td>enable fixed apodization.
%>                                 Apodization window of width
%>                                 @ref n_active_elements and type
%>                                 @ref window with an orientation
%>                                 specified by 3 Euler angles (only
%>                                 supported for bft3_sampled_image)</td></tr>
%>  <tr><td>n_active_elements</td><td>float</td><td>number of active elements (used when fixed = true) </td></tr>
%>  <tr><td>dynamic</td><td>bool</td><td>enable dynamic apodization,
%>                                 specified by an F-number and a
%>                                 window type.</td></tr>
%>  <tr><td>f</td><td>float</td><td>F-number used for dynamic apodization</td></tr>
%>  <tr><td>window</td><td>char</td><td>window function used for
%>                                 dynamic and fixed apodization,
%>                                 'Rectwin', 'Hamming', 'Hann',
%>                                 'Blackman', 'Tukey', 'Gaussian',
%>                                 or 'Bartlett'</td></tr>
%>  <tr><td>window_parameter</td><td>float</td><td>Some windows require a
%>                                 parameter. For the 'Gaussian'
%>                                 window, this is the inverse
%>                                 std. deviation. For the 'Tukey'
%>                                 window, this is ratio of taper to
%>                                 constant sections normalized to
%>                                 (0,1); 0 (Hanning), 1 (Rectwin)</td></tr>
%>  <tr><td>orientation</td><td>float</td><td>Euler angles (not used yet)</td></tr>
%> </table>
%> 
%> @par Read-only properties:
%> <table rules="none">
%>  <tr><td>@ref aperture   </td><td>@ref bft3_aperture</td><td>associated aperture</td></tr>
%>  <tr><td>@ref Handle </td><td>uint(32/64)</td><td>pointer</td></tr>
%>  <tr><td>@ref Id     </td><td>uint(32/64)</td><td>unique identifier</td></tr>
%> </table>
%>
%> @par Methods:
%>    obj @ref bft3_apodization (varargin)\n
%>    obj @ref clone (obj)\n
%>    @ref delete (obj)\n
%>    @ref display (obj)\n
%>
%> $Id: bft3_apodization.m,v 1.59 2012-01-19 10:59:54 jmh Exp $
%>
classdef bft3_apodization < handle

  % Constant properties
  properties (Constant, GetAccess = 'private', SetAccess = ...
              'private')
    %> @brief Library execution unit
    %>
    %> Name of mex-file called by the @ref bft3_apodization class
      
    % Library execution unit
    mexname = 'bft3_mex';%'apodization_mex';
  end

  properties (SetAccess = 'private', GetAccess = 'public')
    %> Object handle
      
    % Object handle
    Handle    
    %> Aperture object
      
    % Aperture object
    aperture
    %> Unique identifier of data content
      
    % Unique identifier of data content
    Id        
  end

  properties (Dependent = true, SetAccess = 'public', GetAccess = 'public')
    %> @brief Apodization reference point
    %>
    %> This point is used as a reference for @ref distances used
    %> for @ref parametric apodization and together with the @ref f
    %> number for the calculation of the size of an active
    %> sub-aperture for @ref dynamic apodization.
      
    % Apodization reference point
    %
    % This point is used as a reference for distances used
    % for parametric apodization and together with the @ref f
    % number for the calculation of the size of an active
    % sub-aperture for dynamic apodization.
    ref               
    %> @brief Parametric apodization enabled
    %>
    %> Enable apodization windows each specified by a distance in
    %> @ref distances and
    %> a window in @ref values. A window is active when we are
    %> further away from @ref ref than the distance.
      
    % Parametric apodization enabled
    %
    % Enable apodization windows each specified by a distance in
    % distances and a window in values. A window is active when we are
    % further away from ref than the distance.
    parametric
    %> @brief Distances to reference point (used when parametric = true)
      
    % Distances to reference point (used when parametric = true)
    distances
    %> @brief float[\#elements, \#dist] Apodization values (used when parametric = true)
    %>
    %> A set of values is specified for each distance used for
    %> parametric apodization. 
      
    % A set of values is specified for each distance used for
    % parametric apodization. 
    values
    %> @brief Dynamic apodization enabled
    %>
    %> This option enables a dynamic apodization using a window
    %> function specified by the @ref window property with a width
    %> calculated using the distance to @ref ref and the F-number @ref f

    % Dynamic apodization enabled
    %
    % This option enables a dynamic apodization using a window
    % function specified by the window property with a width
    % calculated using the distance to ref and the F-number f
    dynamic           
    %> F-number (used when dynamic = true)
      
    % F-number (used when dynamic = true)
    f                 
    %> @brief Window (used when dynamic = true)
    %>
    %> Window function used for dynamic and fixed
    %> apodization. Valid windows are: 'Rectwin', 'Hamming', 'Hann',
    %>                                 'Blackman', 'Tukey', 'Gaussian',
    %>                                 or 'Bartlett'
      
    % Window (used when dynamic = true)
    %
    % Window function used for dynamic and fixed
    % apodization. Valid windows are: 'Rectwin', 'Hamming', 'Hann',
    %                                 'Blackman', 'Tukey', 'Gaussian',
    %                                 or 'Bartlett'
    window            
    %> @brief Window parameter (used when dynamic = true)
    %>
    %> Some windows require a parameter. For the 'Gaussian'
    %>                                 window, this is the inverse
    %>                                 std. deviation. For the 'Tukey'
    %>                                 window, this is ratio of taper to
    %>                                 constant sections normalized to
    %>                                 (0,1); 0 (Hanning), 1 (Rectwin)

    % Window parameter (used when dynamic = true)
    %
    % Some windows require a parameter. For the 'Gaussian'
    %                                 window, this is the inverse
    %                                 std. deviation. For the 'Tukey'
    %                                 window, this is ratio of taper to
    %                                 constant sections normalized to
    %                                 (0,1); 0 (Hanning), 1 (Rectwin)
    window_parameter  
    %> @brief Enable fixed apodization (bft3_sampled image only)
    %>
    %> Apodization window of width @ref n_active_elements and
    %> type @ref window with an @ref orientation specified by 3
    %> Euler angles on the respective aperture
      
    % Enable fixed apodization (bft3_sampled image only)
    %
    % Apodization window of width n_active_elements and
    % type window with an orientation specified by 3
    % Euler angles on the respective aperture 
    fixed             
    %> @brief Number of active elements (used when fixed = true)
    %>
    %> Width of apodization window (when fixed = true), the pitch is taken
    %> as the distance between the first two elements on the
    %> receive aperture, hence it only works for @ref
    %> bft3_aperture's constructed as 'linear array''s
      
    % Number of active elements (used when fixed = true)
    %
    % Width of apodization window (when fixed = true), the pitch is taken
    % as the distance between the first two elements on the
    % receive aperture, hence it only works for 
    % bft3_aperture's constructed as 'linear array''s      
    n_active_elements 
    %> Orientation (used when fixed = true or when the @ref bft3_aperture::ppwave
    %> property is enabled on the aperture)

    % Orientation (used when fixed = true or when the bft3_aperture::ppwave
    % property is enabled on the aperture)
    orientation       
  end

  methods
    % ======================================================================
    %> @brief Class constructor
    %>
    %> Create an apodization object using four arguments or by specifying
    %> a number of options as string-arguments pairs.
    %>
    %> @par Calling:
    %> \n ob = bf3_apodization(aperture, ref, distances, values)        
    %>
    %>  or
    %> @par
    %>          ob = bft3_apodization(aperture, options)
    %>
    %> @par Parameters:
    %> <table rules="none">
    %>  <tr><td>aperture                         </td><td>class bft3_aperture</td><td></td></tr>
    %>  <tr><td>ref is [1 3]                     </td><td>float</td><td>apodization reference point</td></tr>
    %>  <tr><td>distances is [1, # dist]         </td><td>float</td><td>distances to reference point</td></tr>
    %>  <tr><td>values is [# elements, # dist]   </td><td>float</td><td>                          apodization values</td></tr>
    %> </table>                                                         
    %>      
    %> @par Options:
    %> <table rules="none">
    %>  <tr><td>ref is [1,3]                    </td><td>float </td><td>Apodization reference point           </td></tr>
    %>  <tr><td>distances is [1, # dist]        </td><td>float </td><td>Distances to this reference           </td></tr>
    %>  <tr><td>values is [# elements, # dist]  </td><td>float    </td><td>Apodization values (window functions) </td></tr>
    %> </table>
    %>
    %> @param aperture
    %> @param varargin
    %>
    % ======================================================================
    function obj = bft3_apodization(aperture, varargin)

      % Defaults:
      st.ref = [0 0 0];
      st.distances = 0;
      st.dynamic = false;
      st.parametric  = true;
      st.fixed   = false;
      st.window  = 'Hamming';
      st.window_parameter = 1.0;
      st.f       = 1;
      st.orientation = [pi/2 0 0];
      st.n_active_elements = uint32(64);
      
      if (nargin == 0)
        help(mfilename), error('at least one argument is required');
      elseif (nargin == 1)
        if (ischar(aperture) && ...
            strcmp(aperture,'test'))
          % Unit test
          obj = bft3_apodization_test();
          return
        elseif isstruct(aperture)
          % Construct apodization from struct
          st = aperture;
          st.Handle = st.aperture.Handle;
        elseif (length(aperture)==1)
          if (strcmp(class(aperture),'bft3_aperture'))
            % Default constructor
            st.Handle = aperture.Handle;
            st.values = ones(size(aperture.pos,1),1);
          elseif (strcmp(class(aperture),'uint32') || ...
                  strcmp(class(aperture),'uint64'))
            % Apodization returned by line object, must be pointer type
            obj.Handle = aperture;
            return
          else
            help(mfilename), error('at least one argument is required');
          end
        end
      elseif (nargin == 4)
        % Old constructor
        eval(['obj.Handle=', bft3_apodization.mexname,...
              '(''apodization,ctor,manual'',',...
              'aperture.Handle,'...
              'bft3_conv_float(varargin{1}),' ...
              'bft3_conv_float(varargin{2}),' ...
              'bft3_conv_float(varargin{3}));']);
        return
      else
        st.Handle = aperture.Handle;
        st.values = ones(size(aperture.pos,1),1);
        st = bft3_va_arg(st,varargin);
      end
      
      eval(['obj.Handle=', bft3_apodization.mexname,...
            '(''apodization,ctor,manual'',',...
            'st.Handle,'...
            'bft3_conv_float(st.ref),' ...
            'bft3_conv_float(st.distances),' ...
            'bft3_conv_float(st.values));']);

      obj.dynamic           = st.dynamic;
      obj.parametric            = st.parametric;
      obj.fixed             = st.fixed;
      obj.window            = st.window;
      obj.window_parameter  = st.window_parameter;
      obj.f                 = st.f;
      obj.orientation       = st.orientation;
      obj.n_active_elements = st.n_active_elements;
    end
    % ======================================================================
    %> @brief Display function
    %>
    %> Display the properties and member functions of the class
    %>    
    %> @param obj instance of the bft3_apodization class.
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

    % ======================================================================
    %> @brief Clone apodization
    %>
    %> @param obj instance of the bft3_apodization class.
    %> @retval obj instance of the bft3_apodization class (deep copy)
    % ======================================================================
    function out = clone(obj)
      % Clone apodization
      eval(['Handle=', bft3_apodization.mexname,...
            '(''apodization,get,clone'',obj.Handle);']);
      out = bft3_apodization(Handle);
    end

    % ======================================================================
    %> @brief Apodization reference
    %>
    %> @param obj instance of the bft3_apodization class
    %> @retval out reference position
    % ======================================================================
    function out = get.ref(obj)
      % Apodization reference
      eval(['out=', bft3_apodization.mexname,...
            '(''apodization,get,ref'',obj.Handle);']);
    end
    
    % ======================================================================
    %> @brief Distance to reference point (used when parametric = true)
    %>
    %> @param obj instance of the bft3_apodization class
    %> @retval out float[1 # dist]
    % ======================================================================
    function out = get.distances(obj)
      % Distance to reference point (used when parametric = true)
      eval(['out=', bft3_apodization.mexname,...
            '(''apodization,get,distances'',obj.Handle);']);
    end
    % ======================================================================
    %> @brief Apodization values
    %>
    %> @param obj instance of the bft3_apodization class
    %> @retval out float[# elements # dist]
    % ======================================================================
    function out = get.values(obj)
      % Apodization values
      eval(['out=', bft3_apodization.mexname,...
            '(''apodization,get,values'',obj.Handle);']);
    end
    % ======================================================================
    %> @brief Aperture
    %>
    %> @param obj instance of the bft3_apodization class
    %> @retval out instance of the bft3_aperture class
    % ======================================================================
    function out = get.aperture(obj)
      % Aperture
      eval(['Handle =', bft3_apodization.mexname,...
            '(''apodization,get,aperture'',obj.Handle);']);
      out = bft3_aperture(Handle);
    end
    % ======================================================================
    %> @brief Unique identifier
    %>
    %> @param obj instance of the bft3_apodization class
    %> @retval out uint(32/64)
    % ======================================================================
    function out = get.Id(obj)
      % Unique identifier
      eval(['out=', bft3_apodization.mexname,...
            '(''apodization,get,id'',obj.Handle);']);
    end
    % ======================================================================
    %> @brief Dynamic apodization
    %>
    %> @param obj instance of the bft3_apodization class
    %> @retval out bool
    % ======================================================================
    function out = get.dynamic(obj)
      % Dynamic apodization
      eval(['out=', bft3_apodization.mexname,...
            '(''apodization,get,dynamic'',obj.Handle);']);
    end
    % ======================================================================
    %> @brief F-number
    %>
    %> @param obj instance of the bft3_apodization class
    %> @retval out float
    % ======================================================================
    function out = get.f(obj)
      % F-number
      eval(['out=', bft3_apodization.mexname,...
            '(''apodization,get,f'',obj.Handle);']);
    end
    % ======================================================================
    %> @brief Window
    %>
    %> Window can be either 'Rectwin', 'Hamming', 'Hann',
    %> 'Blackman', 'Tukey', 'Gaussian', or 'Bartlett'
    %>
    %> @param obj instance of the bft3_apodization class
    %> @retval out char
    % ======================================================================
    function out = get.window(obj)
      % Window function
      eval(['out=', bft3_apodization.mexname,...
            '(''apodization,get,window'',obj.Handle);']);
    end
    % ======================================================================
    %> @brief Window parameter
    %>
    %> @param obj instance of the bft3_apodization class
    %> @retval out float
    % ======================================================================
    function out = get.window_parameter(obj)
      % Window parameter
      eval(['out=', bft3_apodization.mexname,...
            '(''apodization,get,window_parameter'',obj.Handle);']);
    end
    % ======================================================================
    %> @brief Parametric apodization
    %>
    %> @param obj instance of the bft3_apodization class.
    %> @retval out bool
    % ======================================================================
    function out = get.parametric(obj)
      % Parametric apodization
      eval(['out=', bft3_apodization.mexname,...
            '(''apodization,get,parametric'',obj.Handle);']);
    end
    
    % ======================================================================
    %> @brief Set reference position used for parametric and dynamic apodization
    %>
    %> @param obj instance of the bft3_apodization class.
    %> @param data float[1 3]
    % ======================================================================
    function set.ref(obj,data)
      % Set reference position used for parametric and dynamic apodization
      eval([bft3_apodization.mexname,...
            '(''apodization,set,ref''',...
            ',obj.Handle,bft3_conv_float(data));']);
    end
    % ======================================================================
    %> @brief Set distances used for parametric apodization
    %>
    %> @param obj instance of the bft3_apodization class.
    %> @param data float[1, \#dist]
    % ======================================================================
    function set.distances(obj,data)
      % Set distances used for parametric apodization
      eval([bft3_apodization.mexname,...
            '(''apodization,set,distances''',...
            ',obj.Handle,bft3_conv_float(data));']);
      %      bft3_warn('Some properties have been reset\n')
    end
    % ======================================================================
    %> @brief Set apodization windows for parametric apodization
    %>
    %> @param obj instance of the bft3_apodization class.
    %> @param data float[\#elements, \#dist]
    % ======================================================================
    function set.values(obj,data)
      % Set apodization windows for parametric apodization
      eval([bft3_apodization.mexname,...
            '(''apodization,set,values''',...
            ',obj.Handle,bft3_conv_float(data));']);
      %      bft3_warn('Some properties have been reset\n')
    end
    % ======================================================================
    %> @brief Enable dynamic apodization
    %>
    %> @param obj instance of the bft3_apodization class.
    %> @param data bool
    % ======================================================================
    function set.dynamic(obj,data)
      % Enable dynamic apodization
      eval([bft3_apodization.mexname,...
            '(''apodization,set,dynamic''',...
            ',obj.Handle,data);']);
    end
    % ======================================================================
    %> @brief Set F-number
    %>
    %> @param obj instance of the bft3_apodization class.
    %> @param data float
    % ======================================================================
    function set.f(obj,data)
      % Set F-number
      eval([bft3_apodization.mexname,...
            '(''apodization,set,f''',...
            ',obj.Handle,bft3_conv_float(data));']);
    end
    % ======================================================================
    %> @brief Set parametric apodization
    %>
    %> @param obj instance of the bft3_apodization class.
    %> @param data bool
    % ======================================================================
    function set.parametric(obj,data)
      % Set parametric apodization
      eval([bft3_apodization.mexname,...
            '(''apodization,set,parametric''',...
            ',obj.Handle,data);']);
    end
    % ======================================================================
    %> @brief Set window (used only when dynamic = true)
    %>
    %> Window can be either 'Rectwin', 'Hamming', 'Hann',
    %> 'Blackman', 'Tukey', 'Gaussian', or 'Bartlett'
    %>
    %> @param obj instance of the bft3_apodization class.
    %> @param data char
    % ======================================================================
    function set.window(obj,data)
      % Set window (used only when dynamic = true)
      eval([bft3_apodization.mexname,...
            '(''apodization,set,window''',...
            ',obj.Handle,data);']);
    end
    % ======================================================================
    %> @brief Set window parameter (used only when dynamic = true)
    %>
    %> @param obj instance of the bft3_apodization class.
    %> @param data float
    % ======================================================================
    function set.window_parameter(obj,data)
      % Set window parameter (used only when dynamic = true)
      eval([bft3_apodization.mexname,...
            '(''apodization,set,window_parameter''',...
            ',obj.Handle,bft3_conv_float(data));']);
    end
    % ======================================================================
    %> @brief Enable fixed apodization
    %>
    %> @param obj instance of the bft3_apodization class.
    %> @param data bool
    % ======================================================================
    function set.fixed(obj,data)
      % Enable fixed apodization
      eval([bft3_apodization.mexname,...
            '(''apodization,set,fixed''',...
            ',obj.Handle,data);']);
    end
    % ======================================================================
    %> @brief Fixed apodization
    %>
    %> @param obj instance of the bft3_apodization class.
    %> @retval out bool
    % ======================================================================
    function out = get.fixed(obj)
      % Fixed apodization
      eval(['out=', bft3_apodization.mexname,...
            '(''apodization,get,fixed'',obj.Handle);']);
    end
    % ======================================================================
    %> @brief Set number of active elements
    %>
    %> This property is only used for @ref bft3_sampled_image and
    %> when fixed = true
    %>
    %> @param obj instance of the bft3_apodization class.
    %> @param data uint32
    % ======================================================================
    function set.n_active_elements(obj,data)
      %  Set number of active elements
      %
      % This property is only used for bft3_sampled_image and
      % when fixed = true
      eval([bft3_apodization.mexname,...
            '(''apodization,set,n_active_elements''',...
            ',obj.Handle,data);']);
    end
    % ======================================================================
    %> @brief Number of active elements
    %>
    %> This property is only used for @ref bft3_sampled_image and
    %> when fixed = true
    %> @param obj instance of the bft3_apodization class.
    %> @retval out uint32
    % ======================================================================
    function out = get.n_active_elements(obj)
      % Number of active elements
      %
      % This property is only used for bft3_sampled_image and
      % when fixed = true
      eval(['out=', bft3_apodization.mexname,...
            '(''apodization,get,n_active_elements'',obj.Handle);']);
    end
    % ======================================================================
    %> @brief Set orientation
    %>
    %> Set Euler angles used for orientation of apodization. This
    %> option is not yet implemented
    %>
    %> @param obj instance of the bft3_apodization class.
    %> @param data float[1 3]
    % ======================================================================
    function set.orientation(obj,data)
      % Set orientation
      %
      % Set Euler angles used for orientation of apodization. This
      % option is not yet implemented
      eval([bft3_apodization.mexname,...
            '(''apodization,set,orientation''',...
            ',obj.Handle,bft3_conv_float(data));']);
    end
    % ======================================================================
    %> @brief Orientation
    %>
    %> Euler angles used for orientation of apodization or
    %> plane wave excitations. This option is only valid for
    %> plane-waves and/or beamformation using an @ref
    %> bft3_sampled_image
    %>
    %> @param obj instance of the bft3_apodization class.
    %> @retval out float[1 3]
    % ======================================================================
    function out = get.orientation(obj)
      % Orientation
      %
      % Euler angles used for orientation of apodization or
      % plane wave excitations. This option is only valid for
      % plane-waves and/or beamformation using an
      % bft3_sampled_image
      eval(['out=', bft3_apodization.mexname,...
            '(''apodization,get,orientation'',obj.Handle);']);
    end
    
    % ======================================================================
    %> @brief Class destructor
    %>
    %> Delete method are called before an object of the class is destroyed 
    %> @param obj instance of the bft3_apodization class.
    % ======================================================================
    function delete(obj)
      % Class destructor
      try
        eval([bft3_apodization.mexname,...
            '(''apodization,dtor'',obj.Handle);']);
      catch me
        bft3_warn(['Garbage collector has already freed the memory\n',...
            me.message]);
      end
    end
  end  % methods
end % class

% ======================================================================
%> @brief Unit test of the bft_apodization class
%>
%> Function included for testing consistency
%>
%> @return obj instance of the bft3_apodization class.
% ======================================================================
function obj = bft3_apodization_test()
  recv_aperture = bft3_aperture('type','linear_array');
  no_elements = size(recv_aperture.pos,1);
  focus_point = zeros(1,3);
  distances = [10.0 30.0]/1000;
  values = [hanning(no_elements) hanning(no_elements)];
  obj_old = bft3_apodization(recv_aperture, bft3_conv_float(focus_point),...
                             bft3_conv_float(distances), ...
                             bft3_conv_float(values));
  obj = bft3_apodization(recv_aperture,...
                         'ref',bft3_conv_float(focus_point),...
                         'distances',bft3_conv_float(distances), ...
                         'values',bft3_conv_float(values),...
                         'dynamic',true);
end
