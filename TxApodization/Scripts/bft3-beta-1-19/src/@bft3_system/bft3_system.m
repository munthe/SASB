% System class
%
% Create a system object containing information about sampling
% frequency @ref fs and the speed of sound @ref c. The object is
% created by specifying a number of options as string-argument
% pairs. This object must be created for any beamformation scenario
%
% ob = bft3_system(options)
%
% Options:
%   'fs'           float, sampling rate
%   'c'            float, speed of sound
%
% Properties:
%   'fs'           float, set or get sampling rate
%   'c'            float, set or get speed of sound
%
% Example:
%   fs = 30e6; c = 1480;
%   globals = bft3_system('fs',fs,'c',c);
%  
% $Id: bft3_system.m,v 1.25 2011-08-02 18:53:51 jmh Exp $

% @file bft3_system.m
% @brief System class
%>
%> See class description
% ======================================================================
%> @brief System class
%>
%> Create a system object containing information about sampling
%> frequency @ref fs and the speed of sound @ref c. The object is
%> created by specifying a number of options as string-argument
%> pairs. This object must be created for any beamformation scenario
%>
%> ob = bft3_system(options)
%>
%> @par Options:
%>   'fs'           float, sampling rate\n
%>   'c'            float, speed of sound\n
%>
%> @par Properties:
%>   'fs'           float, set or get sampling rate\n
%>   'c'            float, set or get speed of sound\n
%>
%> @par Example:
%>@code
%>   fs = 30e6; c = 1480;\n
%>   globals = bft3_system('fs',fs,'c',c);
%>@endcode
%>  
%> $Id: bft3_system.m,v 1.25 2011-08-02 18:53:51 jmh Exp $
% ======================================================================
classdef bft3_system < handle
  % Constant properties
  properties (Constant, GetAccess = 'private', SetAccess = ...
              'private')
    %> @brief Library execution unit
      
    % Library execution unit
    mexname = 'bft3_mex';%'aperture_mex';
  end

  properties (SetAccess = 'private', GetAccess = 'private')
    %> @brief uint(32/64) object handle (read-only)

    % uint(32/64) object handle (read-only)
    Handle
  end

  properties (SetAccess = 'public', GetAccess = 'public')
    %> float Sampling frequency
      
    % float Sampling frequency
    fs
    %> float Speed of sound
      
    % float Speed of sound
    c
  end

  properties (SetAccess = 'private', GetAccess = 'public')
    %> char package version
    version
  end

  methods
    % ======================================================================
    %> @brief Class constructor
    %>
    %> Create a system object by specifying a number of options as
    %> string-argument pairs.
    %>
    %> ob = bft3_system(options)
    %>
    %> @par Options:
    %>   'fs'           float, sampling rate\n
    %>   'c'            float, speed of sound\n
    %>
    %> @par Properties:
    %>   'fs'           float, set or get sampling rate\n
    %>   'c'            float, set or get speed of sound\n
    %>
    %> @par Example:
    %>   fs = 30e6; c = 1480;\n
    %>   globals = bft3_system('fs',fs,'c',c);
    %>  
    %> @param varargin a number of options as string-argument pairs
    %>
    %> @retval obj instance of the bft3_system class.
    % ======================================================================
    function obj = bft3_system(varargin)
      % Construct a system object
      %
      %  Calling:  ah = bft3_system(options);
      %
      %  Options:
      %    'fs'             float
      %    'c'              float
    
      % Defaults
      st.fs = 10e6;
      st.c  = 1000;
    
      % Backup
      backup = st;
      if (nargin == 1)
        if (strcmp(varargin{1},'test'))
          % Unit test of system
          obj = bft3_system_test();
          return;
        end
      else
        st = bft3_va_arg(st,varargin);
      end
      
      % Manual system
      pos = zeros(1,3);
      eval(['obj.Handle=', bft3_system.mexname,...
            '(''aperture,ctor,manual'',bft3_conv_float(pos));']);
      
      % Update system according to arguments
      obj.fs = bft3_conv_float(st.fs);
      obj.c = bft3_conv_float(st.c);
      
      if (obj.fs == backup.fs)
        bft3_warn('You are using default fs: %g\n',obj.fs)
      end
      if (obj.c == backup.c)
        bft3_warn('You are using default c: %g\n',obj.c)
      end
    end
  
    % ======================================================================
    %> @brief Display function
    %>
    %> Display the properties and member functions of the class
    %>
    %> @param obj instance of the bft3_system class.
    % ======================================================================
    function display(obj)
      % Display function
      name = inputname(1);
      fprintf('%s is an object of class %s:\n\n', name, ...
              class(obj))
      disp(obj)
    end
    
    % ======================================================================
    %> @brief Sampling frequency
    %>
    %> @param obj instance of the bft3_system class.
    %> @retval out float sampling frequency
    % ======================================================================
    function out = get.fs(obj)
      % Sampling frequency
      eval(['out=', bft3_system.mexname,...
            '(''aperture,get,fs'',obj.Handle);']);
    end

    % ======================================================================
    %> @brief Get the speed of sound
    %>
    %> @param obj instance of the bft3_system class.
    %> @retval out float speed of sound
    % ======================================================================
    function out = get.c(obj)
      % Speed of sound
      eval(['out=', bft3_system.mexname,...
            '(''aperture,get,c'',obj.Handle);']);
    end

    % ======================================================================
    %> @brief Version
    %>
    %> Get the version string for the toolbox matching the CVS tag
    %> used for compilation
    %>
    %> @param obj instance of the bft3_system class.
    %> @retval out char version string
    % ======================================================================
    function out = get.version(obj)
      % Get configuration date
      eval(['out=', bft3_system.mexname,...
            '(''aperture,get,version'',obj.Handle);']);
    end
    
    % ======================================================================
    %> @brief Set sampling frequency
    %>
    %> @param obj instance of the bft3_system class.
    %> @param data float
    % ======================================================================
    function set.fs(obj,data)
      % Set sampling frequency
      % Sampling frequency
      eval([bft3_system.mexname,...
            '(''aperture,set,fs'',obj.Handle,bft3_conv_float(data));']);
    end

    % ======================================================================
    %> @brief Set speed of sound
    %>
    %> @param obj instance of the bft3_system class.
    %> @param data float
    % ======================================================================
    function set.c(obj,data)
      % Set speed of sound
      eval([bft3_system.mexname,...
            '(''aperture,set,c''',...
            ',obj.Handle,bft3_conv_float(data));']);
    end

    % ======================================================================
    %> @brief Class destructor
    %>
    %> Delete method are called before an object of the class is destroyed 
    %> @param obj instance of the bft3_system class.
    % ======================================================================
    function delete(obj)
      % Class destructor
      try
        eval([bft3_system.mexname,...
            '(''aperture,dtor'',obj.Handle);']);
      catch
        % Throws warning if garbage collector for dll has cleaned the object
        warn 'Garbage collector has already freed the memory\n'
      end
    end
  end  % methods
end % class

% ======================================================================
%> @brief Unit test of the bft3_system class
%>
%> Function included for testing consistency. Throws an error in
%> case of in-consistency
%>
%> @return obj instance of the bft3_system class.
% ======================================================================
function obj = bft3_system_test()
  % Unit test of the bft3_system class
  obj = bft3_system('fs', 100e6, 'c', 1480);
end

