function onda_lib_connect(hostname)
%INTERNAL function for the Onda library 2.

port = '49999';
mex_network('init', hostname, port);
resp=mex_network('recv');
if ~strcmp(resp, 'Connected')
    error('Expected ''Connected'' from Onda, got ''%s''', resp);
end
