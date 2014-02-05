function apo = cfu_vs_apo(points, varargin)
%
% apo = cfu_vs_apo(points, paramter-value pairs)
%
% Input parameters:
%  points:     Nx3 vector = [x y z]                             [m]
%  vs:         1x3 coordinate of the virtual source/focus point [m]
%  xdc_center: 1x3 coordinate of aperture center                [m]
%  min_width:  1x1 Minimum width of the apodization (at focus)  [m]
%  alpha:      1x1 Acceptance angle of the apodization          [deg]
%
% Example:
%   points = [linspace(0,50e-3,100)', zeros(100,1) 40*ones(100,1)*1e-3];
%   apo = sa_vs_apo(points, 'vs', [0 0 30e-3], 'alpha',30, 'xdc_center',[-3 -2 0]/1e3)
%
% 2013-09-15, Mofi, Version 1.0
% 2013-09-19, Mofi, Version 1.1, Renamed to cfu_vs_apo.
%
if nargin == 0
    cfu_vs_apo_test;
    name = mfilename;
    help (name)
    return
end

% Input parameters
st.vs           = [0 0 0];
st.xdc_center   = [0 0 0];
st.alpha        = 45;
st.min_width    = 0;
st.debug        = 0;
st.window_type  = 'Hanning';
st.window_param = 0.8;

% Parse Input Parameters
st = cfu_parse_input_parameters(st, varargin);

num_points = size(points,1);

%calc center line
c_line = st.vs - st.xdc_center;
c_line = c_line./norm(c_line);

% distance from VS to point projected onto center line
% dot product
p_proj = [repmat(c_line(1),num_points,1).*points(:,1) + ...
          repmat(c_line(2),num_points,1).*points(:,2) + ...
          repmat(c_line(3),num_points,1).*points(:,3)];
% scale the unit vector
p_center = repmat(c_line,num_points,1) .* repmat(p_proj,1,3);
% distance between the two points (norm of difference)
vs_dist = sqrt( (repmat(st.vs(:,1),num_points,1)-p_center(:,1)).^2  + ...
                (repmat(st.vs(:,2),num_points,1)-p_center(:,2)).^2  + ...
                (repmat(st.vs(:,3),num_points,1)-p_center(:,3)).^2 );

% max allowed distance from line (width of apodization)
max_dist  = tand(st.alpha) * vs_dist;
% make sure apodization is not too narrow
max_dist(max_dist<st.min_width) = st.min_width;

%% distance from point to center line 
% cross product
b = points;
c = points-repmat(c_line, num_points,1);
a = zeros(num_points, 3);
a(:,1) = b(:,2).*c(:,3) - b(:,3).*c(:,2);
a(:,2) = b(:,3).*c(:,1) - b(:,1).*c(:,3);
a(:,3) = b(:,1).*c(:,2) - b(:,2).*c(:,1);
%norm
line_dist = sqrt(a(:,1).^2 + a(:,2).^2 + a(:,3).^2);

% normalised distance
line_dist_n = line_dist./max_dist;

% get the actual apodization
apo = cfu_window(st.window_type, line_dist_n, st.window_param);
apo(line_dist_n>1) = 0;


% Debug plots
if st.debug
    figure;plot(line_dist, '.-');  title('line dist')
    figure;plot(line_dist_n, '.-');title('Normalized line dist')
    figure;plot(max_dist, '.-');   title('Max dist')
    figure;plot(apo, '.-');        title('APO')
end
end





%% Test Function
function cfu_vs_apo_test

vs         = [10e-3 0e-3 30e-3];
alpha      = 25;
xdc_center = [0 0 0]/1e-3;
min_width  = 3e-3;
window_type= 'Tukey';


x_ar = linspace(-100, 100, 1000)/1e3;
y_ar = 0;
z_ar = linspace(0, 120, 1000)/1e3;
[x y z] = meshgrid(x_ar, y_ar, z_ar);
points = [x(:) y(:) z(:)];
apo = cfu_vs_apo(points, 'vs', vs, ...
                'alpha',      alpha, ...
                'xdc_center', xdc_center, ...
                'min_width',  min_width, ... 
                'window_type', window_type);
apo = reshape(apo, 1000, 1000)';
figure;
imagesc(x_ar*1e3, z_ar*1e3, apo)
xlabel('x-axis [mm]')
ylabel('z-axis [mm]')
hold on
plot(vs(1)*1e3, vs(3)*1e3, 'k*', 'MarkerSize', 10, 'LineWidth',2)
plot(xdc_center(1)*1e3, xdc_center(3)*1e3, 'k*', 'MarkerSize', 10, 'LineWidth',2)
axis equal
axis image

x_ar = linspace(-100, 100, 1000)/1e3;
y_ar = linspace(-100, 100, 1000)/1e3;
z_ar = 80e-3;
[x y z] = meshgrid(x_ar, y_ar, z_ar);
points = [x(:) y(:) z(:)];

apo = cfu_vs_apo(points, ...
                'vs',         vs, ...
                'alpha',      alpha, ...
                'xdc_center', xdc_center, ...
                'min_width',  min_width, ... 
                'window_type', window_type);

apo = reshape(apo, 1000, 1000);
figure;
imagesc(x_ar*1e3, y_ar*1e3, apo)
xlabel('x-axis [mm]')
ylabel('y-axis [mm]')
axis equal
axis image



end
