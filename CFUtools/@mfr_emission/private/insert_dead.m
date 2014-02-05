function out = insert_dead(vec,dead,values)
%
% out = insert_dead(vec,dead,values) 
%
%
% Version 1.0, 2012-05-02, MFR
% Version 1.1, 2012-09-02, Now expands the input vector with the length of 'dead', MFR
% Version 1.2  2012-09-13, Inserted value can now be input as third argument. Array can now be expanded.
% 



out = vec;
[m n] = size(vec); %#ok

if nargin < 3
    values = zeros(length(dead), n);
else
    if size(values,1) ~= length(dead), error('Dim1 of ''values'' does not have the correct size.'); end
    if size(values,2) ~= n           , error('Dim2 of ''values does'' not have the correct size.'); end
end

dead2 = sort(dead);
if sum(dead2 ~= dead)>0, error('The ''dead'' array must be sorted.'); end

for idx=0:length(dead)-1
    if  dead(end-idx) <= size(out,1)        % normal case
        out = [out(1:dead(end-idx)-1, :); values(idx+1,:); out(dead(end-idx):end, :)];
    elseif (dead(end-idx)) == size(out,2)+1 % add one extra element
        out = [out; values(idx+1,:)]; %#ok
    else 
        error('Index in ''dead'' is too large (larger than size of ''vec'').');
    end
end


