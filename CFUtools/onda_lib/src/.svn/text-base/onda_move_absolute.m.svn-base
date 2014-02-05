%
% onda_move_absolute(dist [,axis])
%
% dist:  distance in meters to move in absolute coordinates.
%        Must be a vector with three elements representing
%        [dist_x, dist_y, dist_z]. When a scaler is given, ''axis'' 
%        must also be set.
%
% axis:  Optional to use. Can be a vector or scaler. Choses which
%        axis to move. 0=x, 1=y, 2=z, 3=rot1, 4=rot2.
%
% By MFR, 
% Version 1.0, 2012-04-15, Init version.
% Version 1.1, 2013-04-17, Changed to accept single vector argument.
% Version 1.2, 2013-04-22, Re-added ''axis'', now as an optional argument
%
%


function onda_move_absolute(dist, axis)

if nargin == 1
    % assume a vector for movement in (x,y,z) is given
    axis = [0 1 2]; %x,y,z
end
%make sure all dist-elements are numbers
if ~isnumeric(dist)
    error('dist must contain only numbers.')
end
%make sure all axis-elements are numbers
if ~isnumeric(axis)
    error('dist must contain only numbers.')
end
% make sure dist and axis have same length
if length(dist) ~= length(axis)
    error(['axis and dist must contain the same number of ' ...
           'elements.'])
end


% Send the comnmand to the Onda system
for idx = 1:length(axis)
    if axis(idx) >= 0 && axis(idx) <= 2
        dist = dist*1000;  %m to mm convertion
    else
        error('The axis index must lie in the interval [0 2]')
    end
    
    onda_lib('move_absolute', axis(idx), dist(idx));
end


