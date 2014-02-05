function beam_coord = onda_tank2xdc_coord (tank_coord, xdc_n, xdc_origo)
% beam_coord = onda_tank2beam_coord (tank_coord, xdc_n, xdc_origo)
%
% Maps a Onda tank coordinate to a transducer-array coordinate.
%
% tank_coord: Onda tank coordinate, [3x1] vector.
%      xdc_n: Transducer array normal vector (in tank coordinates).
%  xdc_origo: Origo of the transducer array  (in tank coordinates).
%
% 
% By MOFI
% 2013-12-19, Version 1.0, Init Version.
%

% Define helper variable
z_axis = [0;0;1];

% normalise the normal vector
xdc_n  = xdc_n/norm(xdc_n);

% Get axis-angle representation
r = vrrotvec(z_axis,xdc_n);

% Axis-angle to rotation matrix
R = vrrotvec2mat(r);

% Invert the rotation matrix
R = R';  % transpose equals the inverse of a rotation matrix

% Translate coordinate system
translated = tank_coord- xdc_origo;

% Rotate coordinate system
beam_coord = R*translated;

