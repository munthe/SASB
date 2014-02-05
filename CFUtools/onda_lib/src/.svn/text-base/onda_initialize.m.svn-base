%
%
% onda_initialize([host_addr])
%
% Initializes the connection to the Onda system.
%
% host_addr: string containing the IP-address of the host. Defaults to 10.59.44.100.
% 
%
% By MFR, 
% Version 1.0, 2012-04-15, Init version.
% Version 1.1, 2013-04-15, Added check of whether IP address is valid
%




function onda_initialize(host_addr)

IPv4_pattern = '^([01]?\d\d?|2[0-4]\d|25[0-5])\.([01]?\d\d?|2[0-4]\d|25[0-5])\.([01]?\d\d?|2[0-4]\d|25[0-5])\.([01]?\d\d?|2[0-4]\d|25[0-5])$';


if nargin == 0
    onda_lib('init_connection', '10.59.44.100');
else
    if ~ischar(host_addr)
        error('The host address must be a string');
    else
        is_valid_address = ~isempty(regexp(host_addr, IPv4_pattern))
        if ~is_valid_address
            error('The host address appears not to be a valid IPv4 address');
        end
    end
    onda_lib('init_connection', host_addr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% onda_initialise.m ends here
