%
%
% position = onda_get_position
%
% By MFR, 
% Version 1.0, 2012-04-15, Init version.
% Version 1.1, 2013-04-17, Changed to return position for x,y,z axes, always.
%


function position = onda_get_position


pos_x = onda_lib('get_position', 0);
pos_y = onda_lib('get_position', 1);
pos_z = onda_lib('get_position', 2);

position = [pos_x pos_y pos_z]';
position = position/1000;  %mm to m convertion

