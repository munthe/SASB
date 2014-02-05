%
%
% onda_get_high_limit(axis)
%
% axis:      integer choosing which axis to use. (0,1,2 or 'x','y','z')
%
%
% By MFR, 
% Version 1.0, 2013-04-15, Init version.
% Version 1.1, 2013-04-15, Now also accepts strings to identify axes.
%


function limit = onda_get_high_limit(axis)

if nargin ~= 1
    error('This function requires one input.');
end


if strcmp(axis,'x') || axis == 0
    limit = onda_lib('get_high_limit', 0);
elseif strcmp(axis,'y') || axis == 1
    limit = onda_lib('get_high_limit', 1);
elseif strcmp(axis,'z') || axis == 2
    limit = onda_lib('get_high_limit', 2);
else
    warning('Did not identify the axis.')
end

limit = limit/1000;  %mm to m convertion
