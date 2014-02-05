% ob = bft3_image(aperture, aperture, apodizations, apodizations,
%                 lines, options) 
%
% Create an image object using two apertures, two arrays of apodizations,
% an array of lines and a number of options as 
% string-arguments pairs.
%
% Options:
%
% Properties: (can be modified after construction)
%    interp   Options are
%               'nearest' - nearest neighbour interpolation
%               'linear'  - linear interpolation
%               'cubic'   - the four nearest points are used for
%                           fitting a cubic polynomial
%               'spline'  - natural cubic splines
%               'fir'     - upsampling by a factor of 8, using
%                           a predefined LP FIR filter
%                           of order=48
%
%    nthreads                     uint32, number of execution threads or beamformers
%
% Protected properties:
%    Handle
%
%  Methods:
%    beamform(rf_data, delay, i_xmt)  rf_data is float [n_elements n_samples],
%                                      delay [s] is float
%                                      i_xmt is uint32.
%     - Transmission number is only relevant for dynamic transmit apodization.
%
% $Id: bft3_image.m,v 1.32 2011-07-25 15:55:14 jmh Exp $
%
  
% @file bft3_image.m
% @brief Image class
% ======================================================================
%> @brief Image class
%>
%> ob = bft3_image(aperture, aperture, apodizations, apodizations,
%>                 lines, options) 
%>
%> Create an image object using two apertures, two arrays of apodizations,
%> an array of lines and a number of options as 
%> string-arguments pairs.
%>
%> @par Options:
%>
%> @par Properties: (can be modified after construction)
%> <table rules="none">
%>  <tr><td>interp</td><td>char</td><td>'nearest' </td><td>nearest neighbour interpolation      </td></tr>
%>  <tr><td></td><td></td><td>'linear'  </td><td>linear interpolation                 </td></tr>
%>  <tr><td></td><td></td><td>'cubic'   </td><td>the four nearest points are used for fitting a cubic polynomial</td></tr>
%>  <tr><td></td><td></td><td>'spline'  </td><td>natural cubic splines</td></tr>
%>  <tr><td></td><td></td><td>'fir'  </td><td>upsampling by a factor of 8, using
%>                           a predefined LP FIR filter
%>                           of order=48 </td></tr>
%>  <tr><td>nthreads</td><td>uint32</td><td> </td><td>number of
%> execution threads or beamformers</td></tr>
%> </table>
%>
%> @par Read-only properties:
%> <table rules="none">
%>  <tr><td>@ref Handle </td><td>uint(32/64)</td><td>pointer</td></tr>
%> </table>
%>
%> @par Methods:
%> @ref beamform (rf_data, delay, i_xmt)
%> <table rules="none">
%>  <tr><td></td><td></td><td></td></tr>
%>  <tr><td>rf_data</td><td>float[# elements # samples]</td><td></td></tr>
%>  <tr><td>delay</td><td>float</td>                        <td></td></tr>
%>  <tr><td>i_xmt</td><td>uint32</td>
%> <td>Transmission number i_xmt is only used for dynamic
%> transmit apodization.</td></tr>
%> </table>
%>
%> $Id: bft3_image.m,v 1.32 2011-07-25 15:55:14 jmh Exp $
%>
classdef bft3_image < handle
  % Constant properties
  properties (Constant, GetAccess = 'private', SetAccess = ...
              'private')
    %> @brief Library execution unit
    %>
    %> Name of mex-file called by the @ref bft3_image class
      
    % Library execution unit
    mexname = 'bft3_mex';%'image_mex';
  end

  properties (SetAccess = 'private', GetAccess = 'public')
    %> Handle

    % Handle
    Handle
  end % properties
  properties (Dependent = true, SetAccess = 'public', GetAccess = ...
              'public')
    %> @brief Interpolation type
    %> 
    %> Interpolation can be done using either 'nearest', 'linear', 'cubic', 'spline',
    %> or 'fir interpolation

    % Interpolation type
    interp
    %> @brief Number of execution threads
    %>
    %> Number of threads, the scheduler is starting the first
    %> threads in a cyclic order on the cores available on your
    %> system. If any cores are hyperthreaded, they are selected as
    %> the last cores used before multiple threads are executed on
    %> any core.
      
    % Number of execution threads
    nthreads
  end

  methods
    % ======================================================================
    %> @brief Class constructor
    %>
    %> @param xmt_aperture instance of @ref bft3_aperture class
    %> @param rcv_aperture instance of @ref bft3_aperture class
    %> @param xmt_apodizations @ref bft3_apodization[# lines]
    %> @param rcv_apodizations @ref bft3_apodization[# lines]
    %> @param bft_lines @ref bft3_lines[# lines]
    %> @param varargin
    %>
    %> @return instance of the bft3_image class.
    % ======================================================================
    function obj = bft3_image(xmt_aperture, rcv_aperture, xmt_apodizations, rcv_apodizations,bft_lines, varargin)
      if (nargin == 1 && strcmp(xmt_aperture, 'test'))
        obj = bft3_image_test();
      elseif (nargin >  4)
        % Defaults
        st.nthreads = int32(1);
        st.interp   = 'linear';
        
        [recv_ap_handles{1:size(rcv_apodizations,2)}] = rcv_apodizations.Handle;
        [emit_ap_handles{1:size(xmt_apodizations,2)}] = xmt_apodizations.Handle;
        [line_handles{1:size(bft_lines,2)}]           = bft_lines.Handle;
        
        % Convert cell2mat
        line_handles = cell2mat(line_handles);
        emit_ap_handles = cell2mat(emit_ap_handles);
        recv_ap_handles = cell2mat(recv_ap_handles);
          
        eval(['obj.Handle=', bft3_image.mexname,...
              '(''image,ctor,manual'',',...
              'xmt_aperture.Handle,'...
              'rcv_aperture.Handle,'...
              'emit_ap_handles,'...
              'recv_ap_handles,'...
              'line_handles);']);
        
        if (nargin > 5)
          st = bft3_va_arg(st,varargin);
          obj.interp = st.interp;
          obj.nthreads = obj.nthreads;
        end
      else
        help(mfilename), error('at least five arguments are required');
      end
    end
    
    % ======================================================================
    %> @brief Beamform
    %>
    %> @param obj instance of the bft3_image class.
    %> @param rf_data float[# rf_samples # channels]
    %> @param delays  float time of first sample
    %> @param i_xmt uint32 index specifying origin of emission
    %> (used for dynamic transmit apodization only)
    %> @retval out float[# lines # samples]
    % ======================================================================
    function out = beamform(obj, rf_data, delays, i_xmt)
      % Beamform
      eval(['out=', bft3_image.mexname,...
            '(''image,beamform,slow'',obj.Handle''',...
            ',bft3_conv_float(rf_data),bft3_conv_float(delays)''',...
            ',i_xmt);']);
    end
      
    % ======================================================================
    %> @brief Set interpolation type
    %>
    %> data can be either 'nearest', 'linear', 'cubic', 'spline',
    %> or 'fir
    %>
    %> @param obj instance of the bft3_image class.
    %> @param data char 
    % ======================================================================
    function set.interp(obj,data)
      % Set interpolation type, data can be either 'nearest', 'linear',
      % 'cubic', 'spline', or 'fir
      eval([bft3_image.mexname,...
            '(''image,set,interp''',...
            ',obj.Handle,data);']);
    end
    % ======================================================================
    %> @brief Interpolation type
    %>
    %> Return interpolation used
    %>
    %> @param obj instance of the bft3_image class
    %> @retval out char
    % ======================================================================
    function out = get.interp(obj)
      % Interpolation type
      eval(['out=', bft3_image.mexname,...
            '(''image,get,interp'',obj.Handle);']);
    end
    % ======================================================================
    %> @brief Set number of execution threads
    %>
    %> @param obj instance of the bft3_image class.
    %> @param data input
    % ======================================================================
    function set.nthreads(obj,data)
      % Set number of execution threads
      eval([bft3_image.mexname,...
            '(''image,set,nthreads''',...
            ',obj.Handle,int32(data));']);
    end
    % ======================================================================
    %> @brief Number of execution threads
    %>
    %> @param obj instance of the bft3_image class
    %> @retval out
    % ======================================================================
    function out = get.nthreads(obj)
      % Number of execution threads
      eval(['out=', bft3_image.mexname,...
            '(''image,get,nthreads'',obj.Handle);']);
    end
    
    % ======================================================================
    %> @brief Display function
    %>
    %> @param obj instance of the bft3_image class.
    % ======================================================================
    function display(obj)
      % Display function
      name = inputname(1);
      fprintf('%s is an object of class %s:\n\n', name, ...
              class(obj))
      disp(obj)
      fprintf('  Methods:\n')
      bft3_listfunctions(obj);
      fprintf('     interp\tOptions are\n');
      fprintf('\t\t  ''nearest'' - nearest neighbour interpolation\n');
      fprintf('\t\t  ''linear''  - linear interpolation\n');
      fprintf(['\t\t  ''cubic''   - the four nearest points are used ' ...
               'for\n']);
      fprintf('\t\t              fitting a cubic polynomial\n');
      fprintf('\t\t  ''spline''  - natural cubic splines\n');
      fprintf(['\t\t  ''fir''     - upsampling by a factor of 8, using ' ...
               'a predefined LP FIR filter\n']);
      fprintf('\t\t              of order=48\n\n');
    end
      
    % ======================================================================
    %> @brief Class destructor
    %>
    %> Delete method are called before an object of the class is destroyed 
    %> @param obj instance of the bft3_image class.
    % ======================================================================
    function delete(obj)
      try
        eval([bft3_image.mexname,...
            '(''image,dtor'',obj.Handle);']);
      catch me
        bft3_warn(['Garbage collector has already freed the memory\n',...
            me.message]);
      end
    end
   end  % methods
end % class

% subsref.m subsasgn.m

% ======================================================================
%> @brief Unit test of the bft3_image class
%>
%> Function included for testing consistency. Throws an error in
%> case of in-consistency
%>
%> @return obj instance of the bft3_image class
% ======================================================================
function obj = bft3_image_test()

  xmt_aperture = bft3_aperture('type','linear_array');
  rcv_aperture = bft3_aperture('type','linear_array');
  n_elements = size(xmt_aperture.pos,1);
  xmt_apodization = ...
      bft3_apodization(xmt_aperture,...
                       bft3_conv_float(ones(1,3)),...
                       bft3_conv_float(zeros(1,2)),...
                       bft3_conv_float(ones(n_elements,2)));
  rcv_apodization = ...
      bft3_apodization(rcv_aperture,...
                       bft3_conv_float(ones(1,3)),...
                       bft3_conv_float(zeros(1,2)),...
                       bft3_conv_float(ones(n_elements,2)));

  n_lines = 2;
  lines = [];
  ah = [];
  origin = [0 0 0];
  direction = rand(n_lines,3);
  direction = direction ./ repmat(sqrt(sum(direction'.^2))',[1 3]);
  
  dr = 0.2;
  llength = 1;
  
  lines = bft3_lines('origin',origin,'direction',direction,...
                     'dr',dr,'length',llength);
  
  for i=1:n_lines
    ah = [ah bft3_apodization(rcv_aperture, ones(1,3),...
                              zeros(1,2), ones(n_elements,2))];
  end
  obj = bft3_image(xmt_aperture, rcv_aperture, ah, ...
                    ah,lines,'interp','spline');
end
