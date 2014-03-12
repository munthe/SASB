function [matrix,index] = remove_zeros(matrix)
% Removes padding of zeros from matrix and returns stripped matrix and
% index of contents in the original matrix as [x1,y1;x2,y2].
index = zeros(2);
max1 = max(matrix,[],2);
max2 = max(matrix,[],1);
index(1,1) = find(max1~=0,1,'first');
index(2,1) = find(max1~=0,1,'last');
index(1,2) = find(max2~=0,1,'first');
index(2,2) = find(max2~=0,1,'last');
matrix = matrix(index(1,1):index(2,1),index(1,2):index(2,2));
end