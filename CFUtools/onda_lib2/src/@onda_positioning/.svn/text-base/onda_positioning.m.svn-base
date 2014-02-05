% Class for controlling the Onda tank. Note that the PC must be turned on
% and the Onda software running.
classdef onda_positioning < handle
    
    
methods
        
        % Send a command to the Onda system. All commands return a reply.
        function resp=command(obj, cmd)
            cmd_str=sprintf('%s\n', cmd);
            mex_network('send', cmd_str);
            resp=mex_network('recv');
        end
    
    end
    
end

