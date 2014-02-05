function max_idx = cfu_idx_of_max(volume)
%
% Finds the index of the maximum value.
%   max_idx = idx_of_max(volume);
% 
% Output:
%   max_idx.dim1  -- first  dim
%   max_idx.dim2  -- second dim
%   max_idx.dim3  -- third  dim
%
% The function works on both three, two and one dimensional data.
%
% To find the index of the minimum value, one could use the following:
%   min_idx = idx_of_max(-volume);
%
% By Morten F. Rasmussen. 
% Updated April 2012 (now also supports one- and two-dimensional data).
% 2013-08-07, MFR, Rearranged the output order of dim1, dim2 and dim3. 
% 2013-08-16, MFR, Row vectors are no longer interpreted as 2D Matrices.
% 2013-09-17, MFR, Renamed to cfu_sound_speed_calc.
%

dims  = size(volume);
N_dim = length(dims);
if N_dim == 2
    if dims(1) == 1 % then it is a row vector, not a 2D Matrix
        N_dim =1;
    end
end

        
if N_dim == 3
    [max_val_1 max_1] = max(volume,[],1);
    max_surf = squeeze(max_val_1);
    [max_val_2 max_2] = max(max_surf,[],1);
    max_line = squeeze(max_val_2);
    [max_val temp_dim3] = max(max_line,[],2);  %#ok
    temp_dim2 = max_2(temp_dim3);
    temp_dim1 = max_1(1, temp_dim2, temp_dim3);
    
    max_idx.dim1    = temp_dim1;
    max_idx.dim2    = temp_dim2;
    max_idx.dim3    = temp_dim3;
    max_idx.max_val = max_val;
elseif N_dim ==2
    [max_val_1 max_1] = max(volume,[],1);
    [max_val   max_2] = max(max_val_1); %#ok
    temp_dim2       = max_2;
    max_idx.dim1    = max_1(max_2);
    max_idx.dim2    = temp_dim2;
    max_idx.max_val = max_val;
    
elseif N_dim ==1
    [max_val_1 max_1] = max(volume); %#ok
    max_idx.dim1    = max_1;
    max_idx.max_val = max_val_1;
end
