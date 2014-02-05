%
%
% onda_set_low_limit(axis, limit)
%
% axis:      integer choosing which axis to use. (0,1,2 or 'x','y','z')
% limit:     position limit in meters on the chosen axis.
%
%
% By MFR, 
% Version 1.0, 2013-04-15, Init version.
% Version 1.1, 2013-04-15, Now also accepts strings to identify axes.
%


function onda_set_low_limit(axis, limit)

if nargin ~= 2
    error('This function requires two inputs.');
end

limit = limit*1000;  %m to mm convertion

if strcmp(axis,'x') || axis == 0
    onda_lib('set_low_limit', 0, limit);
elseif strcmp(axis,'y') || axis == 1
    onda_lib('set_low_limit', 1, limit);
elseif strcmp(axis,'z') || axis == 2
    onda_lib('set_low_limit', 2, limit);
else
    warning('Did not identify the axis.')
end

