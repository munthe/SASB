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
% Version 1.2, 2013-12-16, Now uses Onda Lib. 2.


function onda_set_position(position)

%make sure all position-elements are numbers
if ~isnumeric(position)
    error('Position must contain only numbers.');
end
if length(position) ~= 3
    error('''position'' must be a vector with 3 elements.');
end

axis = [0 1 2]; %x,y,z
% Send the comnmand to the Onda system
for idx = 1:length(axis)
    if axis(idx) >= 0 && axis(idx) <= 2
        position(idx) = position(idx)*1000;  %m to mm convertion
    end
    
    cmd_str = sprintf('positioner axis%d position %.8f', idx-1, position(idx));
    resp    = onda_lib_command(cmd_str);
    
    if ~strcmp(resp, 'ok')
        error(sprintf('Expected ''ok'' from Onda, but got ''%s''.', resp));
    end
end

