function resp=onda_lib_command(cmd)
% INTERNAL function for the Onda Library 2

cmd_str = sprintf('%s\n', cmd);
mex_network('send', cmd_str);
resp    = mex_network('recv');
