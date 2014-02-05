function idx = row_col2idx (row_tot, col_tot, row, col)
% 
% idx = row_col2idx (row_tot, col_tot, row, col)
%
% returns the index of the chosen column or row, out of the full matrix.
% 
% By Morten R. Rasmussen, 2012-07-10
% Version 1.1, 2013-06-25, Added input check
%
%
%

if max(row(:)) > row_tot || max(col(:)) > col_tot
    error(['''row'' and ''col'' may not contain elements larger than ''row_tot'' or ' ...
           '''col_tot''.'])
end

idx = zeros(row_tot, col_tot);
idx(row, col) = 1;
idx=transpose(idx);
idx=idx(:)';

