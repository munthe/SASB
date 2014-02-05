function id = cfu_get_node_id
% CFU_GET_NODE_ID
%
% Returns ID of the cluster node. If not run on a cluster node, it returns 200.
% on fcfuX returns  X
% on cfuX  returns  100+X
% on mfpc  returns  200
%
% Example:
%   id = get_cluster_id;
%
% Version 1.0, 2012-02-28, Init Version, By Morten F. Rasmussen
% Version 1.1, 2012-03-13, Renamed to mfr_get_node_id,  By Morten F. Rasmussen
% Version 1.2, 2012-09-17, Renamed to cfu_get_node_id,  By Mofi
%

[d name] = system('hostname');

if strcmp(name(1:3), 'fcf')
    id = sscanf(name(5:end), '%i');
elseif strcmp(name(1:3), 'cfu')
    id = sscanf(name(4:end), '%i');
    id=id+100;
else
    id=200;
end


