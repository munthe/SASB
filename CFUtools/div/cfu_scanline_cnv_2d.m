function img_out = cfu_scanline_cnv_2d (img_in, depth_ar, ang_ar, x_sample_pt, y_sample_pt)

% Scan line convert an image
% img_out = cfu_scanline_cnv_2d (img_in, depth_axis,angle_axis [,x_sample_pt, y_sample_pt])
%
% Input:
%   angel_axis must be in radians
%   x_sample_pt is a 1-D array of where to sample the image onthe cartesian x-axis
%   y_sample_pt is a 1-D array of where to sample the image onthe cartesian y-axis
%
% Output:
%   img_out.val   - the pixel values
%   img_out.x     - the x axis in meters
%   img_out.y     - the y axis in meters
%
% img_in may also be the only argument if it is a struct containing all other function
% arguments as members, i.e.:
%   img_in.val    - the pixel values    
%   img_in.depth_axis
%   img_in.angle_axis
%   img_in.x_sample_pt
%   img_in.y_sample_pt
%
% 2011-08-03, MFR, Init version.
% 2013-07-09, MFR, Added the 1D arrays: x_sample_pt and y_sample_pt
%

if nargin == 1
    depth_ar    = img_in.depth_axis;
    ang_ar      = img_in.angle_axis;
    val         = img_in.val;
    x_sample_pt = img_in.x_sample_pt;
    y_sample_pt = img_in.y_sample_pt;
    clear img_in;
    img_in = val;
end

%% cylendrical to cartesian coords
x_coord.val = depth_ar'*sin(ang_ar);
y_coord.val = depth_ar'*cos(ang_ar);
% limits in cartesian
x_coord.min = min(x_coord.val(:));
x_coord.max = max(x_coord.val(:));
y_coord.min = min(y_coord.val(:));
y_coord.max = max(y_coord.val(:));
% sample points in cartesian coords
N_px_x = 800;
N_px_y = 800;
if nargin < 4 || isempty(x_sample_pt)
    img_out.x = linspace(x_coord.min,x_coord.max,N_px_x);
    img_out.y = linspace(y_coord.min,y_coord.max,N_px_y);
else
    img_out.x = x_sample_pt;
    img_out.y = y_sample_pt
end    

[sample.x_grid sample.y_grid] = meshgrid(img_out.x,img_out.y);
% sample points in cylendrical coords
sample.depth  = sqrt(sample.x_grid.^2 + sample.y_grid.^2);
sample.ang    = asin(sample.x_grid ./ sample.depth);

%% scan line convert
%figure;
interp_meth = 'linear';
%interp_meth = 'spline';
img_out.val = interp2 (ang_ar,depth_ar,img_in, ...
                       sample.ang, sample.depth, interp_meth,-400);

