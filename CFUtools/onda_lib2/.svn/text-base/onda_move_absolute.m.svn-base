% Move to a absolute coordinate.
%
%  onda_move_absolute(position[, immediate, tolerance]);
%
% Inputs:
%  position - 1x3 float, coordinate to move to [m],
%  immediate - 0 or 1, return after moving or after sending the command. If
%              set to 1, there is no checking of hitting the limits.
%              Default 0.
%  tolerance - float, maximum allowed tolerance on movement before
%              generating an error [m]. Default (1/92)*1e-3 m.
%
% Version 1.0, 2013-12-16, Init version.
%


function onda_move_absolute(position, immediate, tolerance)

% Input validation
validateattributes(position, {'numeric'}, {'real','nonempty','finite','size', [1 3]});
if nargin<4,    tolerance=(1/92)*1e-3; end
if nargin<3,    immediate=0; end
validateattributes(immediate, {'numeric'}, {'binary', 'scalar'});
validateattributes(tolerance, {'numeric'}, {'positive', 'scalar'});



for idx=1:3
    cmd_str=sprintf('positioner axis%d moveabs %.6f\n', idx-1, position(idx)*1e3);
    resp=onda_lib_command(cmd_str);
    if ~strcmp(resp(1:4), 'None') 
        error('Expected ''None'' from Onda, got ''%s''.', resp);
    end
end




if immediate==0
    pause(0.15); % makes sure Onda has time to respond
    
    for idx=1:4
        end_pos=onda_get_position();
        if sum(isnan(end_pos)) == 0
            break;
        else
            pause(0.15);
        end
    end
    if idx == 4, error('Could not determine the new position.'); end

    
    
    % Error handling
    move_error       = end_pos - position;
    move_error_flag  = abs(move_error)>tolerance;

    if sum(move_error_flag)>0
        fprintf('Movement error exceeds tolerance of %f [mm]\n', tolerance*1e3);
        for idx=1:3
            if  ~move_error_flag(idx)
                fprintf('Axis %d OK.\n', idx-1);
            else
                fprintf('Axis %d NOT OK:\n', idx-1);
                fprintf(' Input position:   %+.3f mm\n', position(idx)*1e3);
                fprintf(' Reached position: %+.3f mm\n', end_pos(idx)*1e3);
                fprintf(' Error:            %+.3f mm\n', move_error(idx)*1e3);
            end
        end
        error(['Movement error exceeds tolerance. You probably reached a hardware (magnet) ' ...
               'limit.']);
    end
end
