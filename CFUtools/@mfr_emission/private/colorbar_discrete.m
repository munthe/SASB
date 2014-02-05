%
% function that creates a colorbar with discrete values 
% By Morten Rasmussen, 2012-12-18
%
% colorbar_discrete(delta_values, [min_val, max_val])
%

%
% Version 1.0 - init version. 
%

function [cb cm] = colorbar_discrete(dval, min_val, max_val)

cvals = caxis;
if nargin < 3
    min_val = cvals(1);
    max_val = cvals(2);
end


values = min_val:dval:max_val;

cm      = colormap;
cm_len  = size(cm,1);
val_len = length(values);

idx    = 1+round(linspace(0,1,val_len)*(cm_len-1));
cm2 = cm(idx,:);
colormap(cm2);
colorbar;
