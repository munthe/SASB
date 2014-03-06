% Line class for Beam Formation
%
% ob = bft3_line(origin, direction, dr, length) 
%
% Create a line object
%
% in:
%    origin        float(1,3)  Starting point of line
%    direction     float(1,3)  Unit vector           
%    dr            float       Size of increment     
%    length        float       Length                
%
% Methods:
%    obj bft3_line(varargin)
%    delete(obj)
%    display(obj)
%    out pos(obj, varargin)
%
% $Id: bft3_line.m,v 1.37 2012-03-26 09:34:07 jmh Exp $
%
  
% @file bft3_line.m
% @brief Line class
% ======================================================================
%> @brief Line class
%>
%> Line class for Beam Formation
%>
%> @par Calling
%> \n ob = bft3_line(origin, direction, dr, length) 
%>
%> @par Parameters:
%> <table rules="none">
%>  <tr><td>origin    </td><td>float[1,3] </td><td>Starting point of line</td></tr>
%>  <tr><td>direction </td><td>float[1,3] </td><td>Unit vector           </td></tr>
%>  <tr><td>dr        </td><td>float      </td><td>Size of increment     </td></tr>
%>  <tr><td>length    </td><td>float      </td><td>Length                </td></tr>
%> </table>
%>
%> @par Methods:
%>    obj bft3_line(varargin)\n
%>    delete(obj)\n
%>    display(obj)\n
%>    out pos(obj, varargin)\n
%>
%> $Id: bft3_line.m,v 1.37 2012-03-26 09:34:07 jmh Exp $
%>
classdef bft3_line < handle
  % TODO:
  %  - Use n_samples instead of length
  
  % Constant properties
  properties (Constant, GetAccess = 'private', SetAccess = ...
              'private')
    %> @brief Library execution unit
    %>
    %> Name of mex-file called by the @ref bft3_line class
      
    % Library execution unit
    mexname = 'bft3_mex';%'line_mex';
  end

  properties (SetAccess = 'private', GetAccess = 'public')
    %> Handle
      
    % Handle
    Handle
    
    %> Origin of line

    % Origin of line
    origin
    
    %> Direction unit vector

    % Direction unit vector
    direction

    %> Line increment
      
    % Line increment
    dr
  end % properties

  properties (Dependent = true, SetAccess = 'private', GetAccess = 'public')
    %> @brief Transmit apodization (read-only)
    %>
    %> (debug purposes only)
      
    % Transmit apodization (read-only)
    %
    % (debug purposes only)
    xmt_apodization
    %> @brief Receive apodization (read-only)
    %>
    %> (debug purposes only)
      
    % Receive apodization (read-only)
    %
    % (debug purposes only)
    rcv_apodization
  end

  methods
    % ======================================================================
    %> @brief Class constructor
    %>
    %> @par Calling
    %> \n ob = bft3_line(origin, direction, dr, length) 
    %>
    %> @par Parameters:
    %> <table rules="none">
    %>  <tr><td>origin    </td><td>float[1,3] </td></tr>
    %>  <tr><td>direction </td><td>float[1,3] </td></tr>
    %>  <tr><td>dr        </td><td>float      </td></tr>
    %>  <tr><td>length    </td><td>float      </td></tr>
    %> </table>
    %>
    %> @param varargin 4 parameters
    %>
    %> @return instance of the bft3_line class.
    % ======================================================================
    function obj = bft3_line(varargin)
      % Construct an object
      if (nargin == 1 && strcmp(varargin{1}, 'test'))
        obj = bft3_line_test();
      elseif (nargin == 4)
        eval(['obj.Handle=', bft3_line.mexname,...
              '(''line,ctor,manual'',',...
              'bft3_conv_float(varargin{1}),' ...
              'bft3_conv_float(varargin{2}),' ...
              'bft3_conv_float(varargin{3}),' ...
              'bft3_conv_float(varargin{4}));']);
      else
        try
          % Defaults:
          st.origin = [0 0 0];
          st.direction = [0 0 1];
          st.dr        = 1;
          st.length    = 100;
          st = bft3_va_arg(st,varargin);
          eval(['obj.Handle=', bft3_line.mexname,...
                '(''line,ctor,manual'',',...
                'bft3_conv_float(st.origin),' ...
                'bft3_conv_float(st.direction),' ...
                'bft3_conv_float(st.dr),' ...
                'bft3_conv_float(st.length));']);
        catch me
          disp(me.message), help(mfilename), error(mfilename)
        end
      end
    end

    function out = get.origin(obj)
      eval(['out=',bft3_line.mexname,...
            '(''line,get,origin'',',...
            'obj.Handle);']);
    end
    
    function out = get.direction(obj)
      eval(['out=',bft3_line.mexname,...
            '(''line,get,direction'',',...
            'obj.Handle);']);
    end
    
    function out = get.dr(obj)
      eval(['out=',bft3_line.mexname,...
            '(''line,get,dr'',',...
            'obj.Handle);']);
    end

    % ======================================================================
    %> @brief Transmit apodization
    %>
    %> @param obj instance of the bft3_line class.
    %> @retval out instance of the bft3_apodization class
    % ======================================================================
    function out = get.xmt_apodization(obj)
      % Transmit apodization
      eval(['Handle=', bft3_line.mexname,...
            '(''line,get,xmt_apodization'',obj.Handle);']);
      if (~isempty(Handle) && ~(ischar(Handle)))
        out = bft3_apodization(Handle);
      else
        out = Handle;
      end
    end

    % ======================================================================
    %> @brief Receive apodization
    %>
    %> @param obj instance of the bft3_line class.
    %> @retval out instance of the bft3_apodization class
    % ======================================================================
    function out = get.rcv_apodization(obj)
      % Receive apodization
      eval(['Handle=', bft3_line.mexname,...
            '(''line,get,rcv_apodization'',obj.Handle);']);

      if (~isempty(Handle) && ~(ischar(Handle)))
        out = bft3_apodization(Handle);
      else
        out = Handle;
      end
    end

    % ======================================================================
    %> @brief Points on the line
    %>
    %> @param obj instance of the bft3_line class.
    %> @param varargin optional indices, size of point are [\#samples 3]
    %> @retval out float[]
    % ======================================================================
    function out = pos(obj,varargin)
      % Points on the line
      eval(['out=', bft3_line.mexname,...
            '(''line,get,pos'',obj.Handle);']);
      if (nargin == 2)
        out = out(varargin{1});
      elseif (nargin == 3)
        out = out(varargin{1},varargin{2});
      end
    end

    % ======================================================================
    %> @brief Transmit apodization values
    %>
    %> Retrieve transmit apodization values for the emission
    %> originating from position index of the transmit aperture
    %>
    %> @param obj instance of the bft3_line class.
    %> @param index of transmit origin
    %> @retval out float[\#samples]
    % ======================================================================
    function out = xmt_apodization_values(obj,index)
      % Retrieve transmit apodization values for the emission
      % originating from position index of the transmit aperture
      eval(['out=', bft3_line.mexname,...
            '(''line,get,xmt_apodization_values'',obj.Handle, uint32(index));']);
    end

    % ======================================================================
    %> @brief Receive apodization values
    %>
    %> Retrieve receive apodization values
    %>
    %> @param obj instance of the bft3_line class.
    %> @retval out float[\#rcv_channels \#samples]
    % ======================================================================
    function out = rcv_apodization_values(obj)
      % Retrieve receive apodization values
      eval(['out=', bft3_line.mexname,...
            '(''line,get,rcv_apodization_values'',obj.Handle);']);
    end
    
    % ======================================================================
    %> @brief Display function
    %>
    %> Display the properties and member functions of the class
    %>    
    %> @param obj instance of the bft3_line class
    % ======================================================================
    function display(obj)
      % Display the properties and member functions of the class
      name = inputname(1);
      fprintf('%s is an object of class %s:\n\n', name, ...
              class(obj))
      disp(obj)
      fprintf('  Methods:\n')
      bft3_listfunctions(obj);
      fprintf('\n');
    end
    
    % ======================================================================
    %> @brief Class destructor
    %>
    %> Delete method are called before an object of the class is destroyed 
    %> @param obj instance of the bft3_image class.
    % ======================================================================
    function delete(obj)
      % Class destructor
      try
        eval([bft3_line.mexname,...
            '(''line,dtor'',obj.Handle);']);
      catch
        bft3_warn 'Garbage collector has already freed the memory\n'
      end
    end
  end  % methods
end % class

% ======================================================================
%> @brief Unit test of the bft_line class
%>
%> Function included for testing consistency. Throws an error in
%> case of in-consistency
%>
%> @return obj instance of the bft3_line class.
% ======================================================================
function obj = bft3_line_test()
  n_elements = 9;
  n_lines    = 9;
  n_samples  = 9;
  
  xmt_aperture = bft3_aperture('type','linear_array','n_elements',9);
  rcv_aperture = bft3_aperture('type','linear_array','n_elements',9);
  
  nine_lines = bft3_lines('origin',xmt_aperture.pos,...
                          'direction',[0 0 1],...
                          'dr', 1/1000,...
                          'length', 9/1000);

  xmt_apo = zeros(n_lines,n_samples); % 
  for i_line=1:n_lines
    rcv_apodization = bft3_apodization(rcv_aperture,'ref',nine_lines(i_line).pos(1,:));
    rcv_apodization.parametric = false;
    rcv_apodization.dynamic    = true;
    rcv_apodization.f          = 0.5;
    xmt_apodization = bft3_apodization(xmt_aperture,'ref',nine_lines(i_line).pos(1,:));
    xmt_apodization.parametric = false;
    xmt_apodization.dynamic    = true;
    xmt_apodization.f          = 0.5;
    bft3_mex('line,set,xmt_apodization',nine_lines(i_line).Handle,...
             xmt_apodization.Handle);
    bft3_mex('line,set,rcv_apodization',nine_lines(i_line).Handle,...
             rcv_apodization.Handle);
    xmt_apo(i_line,:) = nine_lines(i_line).xmt_apodization_values(5);
  end

  imagesc(xmt_apo)
  figure;
  imagesc(nine_lines(5).rcv_apodization_values)
  
  % Receive and transmit apodization for this setup must satisfy
  if (max(abs(nine_lines(5).rcv_apodization_values - xmt_apo)) > ...
      eps)
    error('Error in apodization calculation');
  end
  obj = nine_lines(5);
end
