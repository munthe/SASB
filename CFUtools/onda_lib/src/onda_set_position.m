%
%
% onda_set_position(position)
%
% position:  Vector with three elements setting the position of the
%            x-, y- and z-axis. In units of meters.
%
%
% By MFR, 2012-04-15
% Version 1.0 Init version.
% Version 1.1, 2013-04-17, Changed to set the position for x,y,z axes, always.
%


function onda_set_position(position)

%make sure all position-elements are numbers
if ~isnumeric(position)
    error('dist must contain only numbers.');
end
if length(position) ~= 3
    error('''position'' must be a vector with 3 elements.');
end

axis = [0 1 2]; %x,y,z
% Send the comnmand to the Onda system
for idx = 1:length(axis)
    if axis(idx) >= 0 && axis(idx) <= 2
        position = position*1000;  %m to mm convertion
    end
    
    onda_lib('set_position', axis(idx), position(idx));
end



