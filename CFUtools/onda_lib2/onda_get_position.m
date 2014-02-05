% Get the current position in the Onda coordinate system. This system can
% be defined with set_position.
%
%  position=my_onda.get_position();
%
% Outputs:
%  position - 1x3 float, the current position in [m]

function position=get_position(obj)
    position=zeros(1, 3);
    for idx=1:3
        cmd_str = sprintf('positioner axis%d position?', idx-1);
        resp    = onda_lib_command(cmd_str);
        %TODO: make sure 'resp' is a number.
        position(idx)=str2double(resp)*1e-3;
    end
end

