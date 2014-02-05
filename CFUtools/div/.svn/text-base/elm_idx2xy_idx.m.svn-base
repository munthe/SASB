function [xidx yidx] = elm_idx2xy_idx (elm_idx)
%
% [xidx yidx] = elm_idx2xy_idx (elm_idx)
%
% Version 0.1, Init version, 2012-10-03, mofi.
% Version 0.2, Changed function call from vol_idx_of max to idx_of_max, 2013-06-20, mofi.

% a very ugly hack..  should just be calculated.
T = zeros(32,32);
T(elm_idx)=1;
dim = cfu_idx_of_max(T);
xidx = dim.dim1;
yidx = dim.dim2;



