function tank_coord = onda_xdc2tank_coord (xdc_coord, xdc_n, xdc_origo)
% tank_coord = onda_xdc2tank_coord (xdc_coord, xdc_n, xdc_origo)
%
% Maps a xdc coordinate to a Onda tank coordinate.
%
% xdc_coord: Coordinate to be mapped from the transducer-array coordinate system to the
%            Onda tank coordinate system.    It must be a [3x1] vector.
%     xdc_n: Transducer array normal vector. It must be a [3x1] vector. 
% xdc_origo: Origo of the transducer array.  It must be a [3x1] vector.
%
% xdc_n and xdc_origo must be in Onda tank coordinates.
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

% rotate coordinate system
rotated = R*xdc_coord;


% Translate coordinate system
tank_coord = rotated + xdc_origo;

