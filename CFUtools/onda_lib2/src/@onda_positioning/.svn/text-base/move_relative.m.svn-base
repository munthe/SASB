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

function move_relative(obj, amount, immediate, tolerance)
    validateattributes(amount, {'numeric'}, {{'size', [1 3]}});
    if nargin<4
        tolerance=1/92000;
    end
    if nargin<3
        immediate=0;
    end
    validateattributes(immediate, {'numeric'}, {'binary', 'scalar'});
    validateattributes(tolerance, {'numeric'}, {'positive', 'scalar'});
    
    if immediate==0
        start_pos=obj.get_position();
    else
        warning('Non-immediate mode not implemented');
        % Note that for immediate mode, the resulting 'None' is not flushed
        % from the receive buffer leading to asynchrony between commands
        % and responses.
    end
    
    for idx=1:3
        cmd_str=sprintf('positioner axis%d moverel %f\n', idx-1, amount(idx));
        resp=obj.command(cmd_str);
        if ~strcmp(resp, 'None')
            error('Expected ''None'' from Onda, got %s', resp);
        end
    end
    
    if immediate==0
        end_pos=obj.get_position();
    end
    
    move_error=abs(end_pos-start_pos);
    move_ok=move_error>tolerance;
    
    if sum(move_ok)>0
        fprintf('Movement error exceeds tolerance of %f m', tolerance);
        for idx=1:3
            if move_ok(idx)
                fprintf('Axis %d OK\n', idx-1);
            else
                fprintf('Axis %d NOT OK\n', idx-1);
                fprintf(' Start position: %f\n', start_pos(idx));
                fprintf(' End position:   %f\n', end_pos(idx));
                fprintf(' Error:          %f\n', move_error(idx));
            end
        end
        error('Movement error exceeds tolerance');
    end
end
