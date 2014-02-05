% Move a relative amount.
%
%  my_onda.move_relative(amount[, immediate, tolerance]);
%
% Inputs:
%  amount - 1x3 float, amount to move [m],
%  immediate - 0 or 1, return after moving or after sending the command. If
%              set to 1, there is no checking of hitting the limits.
%              Default 0.
%  tolerance - float, maximum allowed tolerance on movement before
%              generating an error [m]. Default 1/92000.
%
% Version 1.1, 2013-12-16, Now uses procedural calls instead of OOP.
%


function onda_move_relative(amount, immediate, tolerance)

% Input validation
validateattributes(amount, {'numeric'}, {'real','nonempty','finite','size', [1 3]});
if nargin<4,    tolerance=1/92000; end
if nargin<3,    immediate=0; end
validateattributes(immediate, {'numeric'}, {'binary', 'scalar'});
validateattributes(tolerance, {'numeric'}, {'positive', 'scalar'});


if immediate==0
    start_pos=onda_get_position();
end

for idx=1:3
    cmd_str=sprintf('positioner axis%d moverel %f\n', idx-1, amount(idx)*1e3);
    resp=onda_lib_command(cmd_str);
    if ~strcmp(resp, 'None')
        error('Expected ''None'' from Onda, got %s', resp);
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
    positioner_moved = end_pos-start_pos;
    move_error       = positioner_moved - amount;
    move_error_flag  = abs(move_error)>tolerance;

    if sum(move_error_flag)>0
        fprintf('Movement error exceeds tolerance of %f [mm]\n', tolerance*1e3);
        for idx=1:3
            if  ~move_error_flag(idx)
                fprintf('Axis %d OK\n', idx-1);
            else
                fprintf('Axis %d NOT OK\n', idx-1);
                fprintf(' Start position: %.3f mm\n', start_pos(idx)*1e3);
                fprintf(' End position:   %.3f mm\n', end_pos(idx)*1e3);
                fprintf(' Error:          %.3f mm\n', move_error(idx)*1e3);
            end
        end
        error('Movement error exceeds tolerance. ');
    end
end
