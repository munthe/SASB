function pid_out = cfu_log_pid
%
% pid = cfu_log_pid
%
% This function writes MATLAB's PID (Process ID) to the STD ERR file stream.
%
%
% 2013-09-16, MOFI, Init version.
%

% Get the process ID
pid   = feature('getpid');

% Print the PID
stderr = 2;
fprintf(stderr, 'PID: %i\n', pid);

% set the output
if nargout > 0
    pid_out = pid;
end