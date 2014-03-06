% ob = bft3_sampled_image(aperture, aperture, im_geom)
%
% Create an image object using two apertures, an im_geom object,
% and a number of options as string-arguments pairs.
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
%    nthreads                     uint32, number of threads or beamformers
%
% Protected properties:
%    Handle
%    Id
%
%  Methods:
%
% $Id: bft3_sampled_image.m,v 1.20 2011-04-27 20:35:28 jmh Exp $
%

% @file bft3_sampled_image.m
% @brief Sampled Image class
% ======================================================================
%> @brief Sampled Image class
%>
%> ob = bft3_sampled_image(aperture, aperture, im_geom)
%>
%> Create an image object using two apertures, an im_geom object,
%> and a number of options as string-arguments pairs.
%>
%> Options:
%>
%> Properties: (can be modified after construction)
%>    interp   Options are
%>               'nearest' - nearest neighbour interpolation
%>               'linear'  - linear interpolation
%>               'cubic'   - the four nearest points are used for
%>                           fitting a cubic polynomial
%>               'spline'  - natural cubic splines
%>               'fir'     - upsampling by a factor of 8, using
%>                           a predefined LP FIR filter
%>                           of order=48
%>
%>    nthreads                     uint32, number of threads or beamformers
%>
%> Protected properties:
%>    Handle
%>    Id
%>
%>  Methods:
%>
%> $Id: bft3_sampled_image.m,v 1.20 2011-04-27 20:35:28 jmh Exp $
%>
classdef bft3_sampled_image < handle
  % Constant properties
  properties (Constant, GetAccess = 'private', SetAccess = ...
              'private')
    %> @brief Library execution unit
    %>
    %> Name of mex-file called by the @ref bft3_sampled_image class
      
    % Library execution unit
    mexname = 'bft3_mex';%'sampled_image_mex';
  end

  properties (SetAccess = 'private', GetAccess = 'public')
    %> @brief Handle

    % Handle
    Handle
  end
  properties (Dependent = true, SetAccess = 'public', GetAccess = 'public')
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
    %> @param varargin a number of options as string-argument pairs
    %>
    %> @return instance of the bft3_sampled_image class.
    % ======================================================================
    function obj = bft3_sampled_image(varargin)
      if (nargin == 1 && strcmp(varargin{1}, 'test'))
        obj = bft3_sampled_image_test();
      elseif (nargin ~= 5)
          help(mfilename), error(mfilename)
      else
        eval(['obj.Handle=',bft3_sampled_image.mexname,...
              '(''sampled_image,ctor'',',...
              'varargin{1}.Handle,'...
              'varargin{2}.Handle,'...
              'varargin{3}.Handle,'...
              'varargin{4}.Handle, varargin{5});']);
      end
    end
    
    % ======================================================================
    %> @brief Beamform
    %>
    %> @param obj instance of the bft3_sampled_image class.
    %> @param rf_data
    %> @param delay
    %> @param i_xmt
    %> @param angles
    %>
    %> @retval 
    % ======================================================================
    function out = beamform(obj, rf_data, delay,...
                            i_xmt, angles)
      eval(['out=',bft3_sampled_image.mexname,...
            '(''sampled_image,beamform,slow'',',...
            'obj.Handle, bft3_conv_float(rf_data),'...
            'bft3_conv_float(delay),i_xmt,'...
            'bft3_conv_float(angles));']);
    end
      
    % ======================================================================
    %> @brief Set interpolation type
    %>
    %> @param obj instance of the bft3_sampled_image class.
    %> @param data input
    % ======================================================================
    function set.interp(obj,data)
      % Set interpolation type, data can be either 'nearest', 'linear',
      % 'cubic', 'spline', or 'fir
      eval([bft3_sampled_image.mexname,...
            '(''sampled_image,set,interp''',...
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
      eval(['out=', bft3_sampled_image.mexname,...
            '(''sampled_image,get,interp'',obj.Handle);']);
    end
    % ======================================================================
    %> @brief Set number of execution threads
    %>
    %> @param obj instance of the bft3_sampled_image class.
    %> @param data input
    % ======================================================================
    function set.nthreads(obj,data)
      % Set number of execution threads
      eval([bft3_sampled_image.mexname,...
            '(''sampled_image,set,nthreads''',...
            ',obj.Handle,int32(data));']);
    end
    % ======================================================================
    %> @brief Get number of execution threads
    %>
    %> @param obj instance of the bft3_sampled_image class
    %> @retval out
    % ======================================================================
    function out = get.nthreads(obj)
      % Number of execution threads
      eval(['out=', bft3_sampled_image.mexname,...
            '(''sampled_image,get,nthreads'',obj.Handle);']);
    end
    % ======================================================================
    %> @brief Display function
    %>
    %> @param obj instance of the bft3_sampled_image class.
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
      % fprintf(['     interp\tOptions are\n']);
      % fprintf(['\t\t  ''nearest'' - nearest neighbour interpolation\n']);
      % fprintf(['\t\t  ''linear''  - linear interpolation\n']);
      % fprintf(['\t\t  ''cubic''   - the four nearest points are used ' ...
      %          'for\n']);
      % fprintf(['\t\t              fitting a cubic polynomial\n']);
      % fprintf(['\t\t  ''spline''  - natural cubic splines\n']);
      % fprintf(['\t\t  ''fir''     - upsampling by a factor of 8, using ' ...
      %          'a predefined LP FIR filter\n']);
      % fprintf(['\t\t              of order=48\n\n']);
    end
      
    % ======================================================================
    %> @brief Class destructor
    %>
    %> Delete method are called before an object of the class is destroyed 
    %> @param obj instance of the bft3_sampled_image class.
    % ======================================================================
    function delete(obj)
      try
        eval([bft3_sampled_image.mexname,...
            '(''sampled_image,dtor'',obj.Handle);']);
      catch
        bft3_warn 'Garbage collector has already freed the memory\n'
      end
    end
   end  % methods
end % class

% subsref.m subsasgn.m

% ======================================================================
%> @brief Test of bft3_sampled_image class
%>
%> Function included for testing consistency
%>
%> @return obj instance of the bft3_sampled_image class
% ======================================================================
function obj = bft3_sampled_image_test()
  im = bft3_im_geom('test');
  xmt_aperture    = bft3_aperture('test');
  rcv_aperture    = bft3_aperture('test');
  xmt_apodization = bft3_apodization(xmt_aperture,[0 0 0], 0, ones(64,1));
  rcv_apodization = bft3_apodization(rcv_aperture,[0 0 0], 0, ones(64,1));
  
  obj = bft3_sampled_image(xmt_aperture, rcv_aperture,...
                           xmt_apodization, rcv_apodization,im);
end
