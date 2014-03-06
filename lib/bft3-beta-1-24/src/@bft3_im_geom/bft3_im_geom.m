% Image Geometry Class
%
% Create a "image geometry" class that describes the sampling
% characteristics of a single 2D image.
%
% options:
%	'nx'		image dimension
%	'nz'		image dimension (default: nx)
%	'dx'		pixel size (required)
%	'dz'		pixel size (default: -dx)
%	'offset_x'	[units of dx] (default: 0)
%	'offset_z'	[units of dz] (default: 0)
%	'fov'		nx * dx
%	'mask'		logical support mask
%
% out:
%
% methods:
%	x	1D x coordinates of each pixel (1D)
%	z	1D z ""
%	np	sum(mask(:)) (# of pixels to be estimated)
%    circ    2D image with ellipsis
%
%              |          (3)
%    ----------o---------------|-------------------------->
%              |  dx     offset_x
%              |
%              |
%              |
%              |
%              - offset_y
%              |           
%              |           
%              |
%              v
%
%
% dx and dy specifies direction of axes
%
% $Id: bft3_im_geom.m,v 1.13 2011-08-04 18:18:04 jmh Exp $

% @file bft3_im_geom.m
% ======================================================================
%> @brief Image Geometry Class
%>
%> Create a "image geometry" class that describes the sampling
%> characteristics of a single 2D image.
%>
%> @par Options:
%> <table rules="none">
%>  <tr><td>@b 'nx'</td><td>float</td><td>image dimension               </td></tr>
%>  <tr><td>@b 'nz'</td><td>float</td><td>image dimension (default: nx) </td></tr>
%>  <tr><td>@b 'dx'</td><td>float</td><td>pixel size (required)         </td></tr>
%>  <tr><td>@b 'dz'</td><td>float</td><td>pixel size (default: -dx)     </td></tr>
%>  <tr><td>@b 'offset_x'</td><td>float</td><td>[units of dx] (default: 0)</td></tr>
%>  <tr><td>@b 'offset_z'</td><td>float</td><td>[units of dz] (default: 0)</td></tr>
%>  <tr><td>@b 'fov'</td><td>float</td><td>nx * dx </td></tr>
%>  <tr><td>@b 'mask'</td><td>float</td><td>logical support mask </td></tr>
%> </table>
%>
%> @par Methods:
%> <table rules="none">
%>  <tr><td>x</td><td> float[1 np]</td><td>1D x coordinates of each pixel</td></tr>
%>  <tr><td>y</td><td> float[1 np]</td><td>1D y </td></tr>
%>  <tr><td>z</td><td> float[1 np]</td><td>1D z </td></tr>
%>  <tr><td>np</td><td> float</td><td>sum(mask(:)) (\# of pixels to be
%> estimated) </td></tr>
%>  <tr><td>circ </td><td> float[nx nz] </td><td>2D image with ellipsis</td></tr>
%> </table>
%>
%>@verbatim
%>              |          (3)
%>    ----------o---------------|-------------------------->
%>              |  dx     offset_x
%>              |
%>              |
%>              |
%>              |
%>              - offset_z
%>              |           
%>              |           
%>              |
%>              v
%>
%> dx and dz specifies direction of axes
%>@endverbatim
%>
%> $Id: bft3_im_geom.m,v 1.13 2011-08-04 18:18:04 jmh Exp $
%>
classdef bft3_im_geom  < handle

  properties (SetAccess = 'private', GetAccess = 'private')
  end
  properties (SetAccess = 'public', GetAccess = 'public')
    %> Number of x-pixels
    nx        
    %> Number of y-pixels
    ny
    %> Number of z-pixels
    nz        
    %> Pixel separation, x-direction
    dx        
    %> Pixel separation, y-direction
    dy        
    %> Pixel separation, z-direction
    dz        
    %> Pixel offset in units of dx
    offset_x  
    %> Pixel offset in units of dy
    offset_y  
    %> Pixel offset in units of dz
    offset_z  
    %> Field of view (FOV)
    fov       
    %> Image mask
    mask      
    %> Dimensions
    dim       
  end % properties

  methods
    % ======================================================================
    %> @brief Class constructor
    %>
    %>
    %> Create an "image geometry" object that describes the sampling
    %> characteristics of a single 2D image. The object is
    %  constructed by specifying a number of options as 
    %> string-arguments pairs.
    %>
    %> @par Options:
    %> <table rules="none">
    %>  <tr><td>@b 'nx'</td><td>float</td><td>image dimension               </td></tr>
    %>  <tr><td>@b 'nz'</td><td>float</td><td>image dimension (default: nx) </td></tr>
    %>  <tr><td>@b 'dx'</td><td>float</td><td>pixel size (required)         </td></tr>
    %>  <tr><td>@b 'dz'</td><td>float</td><td>pixel size (default: -dx)     </td></tr>
    %>  <tr><td>@b 'offset_x'</td><td>float</td><td>[units of dx] (default: 0)</td></tr>
    %>  <tr><td>@b 'offset_z'</td><td>float</td><td>[units of dz] (default: 0)</td></tr>
    %>  <tr><td>@b 'fov'</td><td>float</td><td>nx * dx </td></tr>
    %>  <tr><td>@b 'mask'</td><td>float</td><td>logical support mask </td></tr>
    %> </table>
    %>
    %> @par varargin
    % ======================================================================
    function obj = bft3_im_geom(varargin)
      % Construct an object    
      if nargin == 1 && strcmp(varargin{1}, 'test')
        obj = bft3_im_geom_test; 
        return;
      end
      if nargin < 1, help(mfilename), error(mfilename), end
      
      % defaults
      obj.nx = [];
      obj.ny = [];
      obj.nz = [];
      obj.dx = [];
      obj.dy = [];
      obj.dz = [];
      obj.offset_x = 0;
      obj.offset_y = 0;
      obj.offset_z = 0;
      obj.fov = [];
      obj.mask = [];
      
      obj = bft3_va_arg(obj, varargin);
      
      % dimensions
      if isempty(obj.ny)
        if isempty(obj.nz)
          obj.ny = obj.nx;
        end
      end
        
      % distances
      if isempty(obj.dx)
        if isempty(obj.fov)
          error 'dx or fov required'
        end
        obj.dx = obj.fov / obj.nx;
      elseif isempty(obj.fov) && ~isempty(obj.nx)
        obj.fov = obj.nx * obj.dx;
      end
      if obj.fov ~= obj.nx * obj.dx
        error 'bad fov'
      end
      
      if isempty(obj.ny)
        obj.dim = [obj.nz obj.nx];
      else
          obj.dim = [obj.nx obj.ny];
      end
      
      if ~isempty(obj.ny)
        if isempty(obj.dy)
          obj.dy = -obj.dx;
        end
      end
      
      % mask
      if isempty(obj.mask)
        obj.mask = true(obj.dim);
      elseif ~islogical(obj.mask)
        error 'mask must be logical'
      elseif isempty(obj.nz)
        if ndims(obj.mask) ~= length(obj.dim)  ...
              || size(obj.mask,1) ~= obj.nx ...
              || size(obj.mask,2) ~= obj.ny
          size(obj.mask), obj.nx, obj.ny
          error 'bad input mask size'
        end
      elseif ndims(obj.mask) ~= length(obj.dim)  ...
            || size(obj.mask,1) ~= obj.nx ...
            || size(obj.mask,2) ~= obj.nz
        size(obj.mask), obj.nx, obj.nz
        error 'bad input mask size'
      end
    end
    
    % ======================================================================
    %> @brief Display function
    %>
    %> @param obj instance of the bft3_im_geom class
    % ======================================================================
    function display(obj)
      name = inputname(1);
      fprintf('%s is an object of class %s:\n\n', name, ...
              class(obj))
      disp(obj)
      fprintf('  Methods:\n')
      bft3_listfunctions(obj);
      fprintf('\n');
        % fprintf('\tx:\t\t return x-coordinates\n')
        % fprintf('\ty:\t\t return y-coordinates\n')
        % fprintf('\tnp:\t\t return number of pixels inside mask\n')
        % fprintf('\tcirc(r):\t return circular image with radius r\n\n')
    end
    % ======================================================================
    %> @brief x coordinates
    %>
    %> @param obj instance of the bft3_im_geom class
    %> @param varargin indices (if any)
    % ======================================================================
    function x = x(obj,varargin)
        % x-coordinates
        wx = (obj.nx-1)/2 + obj.offset_x;
        x = ((0:obj.nx-1)' - wx) * obj.dx;
        x = x(varargin{:});
    end
        

    % ======================================================================
    %> @brief y coordinates
    %>
    %> @param obj instance of the bft3_im_geom class
    %> @param varargin indices (if any)
    % ======================================================================
    function y = y(obj, varargin)
        % y-coordinates
        wy = (obj.ny-1)/2 + obj.offset_y;
        y = ((0:obj.ny-1)' - wy) * obj.dy;
        y = y(varargin{:});
    end

    % ======================================================================
    %> @brief z coordinates
    %>
    %> @param obj instance of the bft3_im_geom class
    %> @param varargin indices (if any)
    % ======================================================================
    function z = z(obj, varargin)
        % z-coordinates
        wz = (obj.nz-1)/2 + obj.offset_z;
        z = ((0:obj.nz-1)' - wz) * obj.dz;
        z = z(varargin{:});
    end
    
    % ======================================================================
    %> @brief Number of pixels
    %>
    %> @param obj instance of the bft3_im_geom class
    % ======================================================================
    function np = np(obj)
      % Number of pixels inside mask
      np = sum(obj.mask(:));
    end

    % ======================================================================
    %> @brief Circle for masking
    %>
    %> @param obj instance of the bft3_im_geom class
    % ======================================================================
    function circ = circ(obj, rad, over)
      % Circle centered at (0,0)
      if ~isvar('rad') || isempty(rad)
	rad = min(abs((obj.nx/2-1)*obj.dx), abs((obj.ny/2-1)*obj.dy));
      end
      if ~isvar('over') || isempty(over)
        over=4;
      end
      
      theta = 0; cx = 0; cy = 0;
      circ = zeros(obj.nx*over,obj.ny*over);
      wx = (obj.nx*over-1)/2 + obj.offset_x*over;
      wy = (obj.ny*over-1)/2 + obj.offset_y*over;
      xx = ((0:obj.nx*over-1) - wx) * obj.dx / over;
      yy = ((0:obj.ny*over-1) - wy) * obj.dy / over;
      [xx yy] = ndgrid(xx, yy);
      x = cos(theta) * (xx-cx) + sin(theta) * (yy-cy);
      y = -sin(theta) * (xx-cx) + cos(theta) * (yy-cy);
      tmp = (x / rad).^2 + (y / rad).^2 <= 1;
      circ(tmp>0) = 1;
      circ = imresize(circ,[obj.nx obj.ny],'bilinear');
    end
  end % methods
end

% ======================================================================
%> @brief Unit test of the bft_im_geom class
%>
%> Function included for testing consistency
%>
%> @return obj instance of the bft3_im_geom class.
% ======================================================================
function obj = bft3_im_geom_test
  %
  % bft3_im_geom_test()
  %
  obj = bft3_im_geom('dx', 2, 'nx', 8, 'dz', 1, 'nz', 8, 'offset_z', -4);
end
