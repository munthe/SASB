% Function to quantize the apodization vector
%
% Inputs:
% Matrix to be quantized
% Quantization levels
%
% Outputs:
% Quantized matrix

function [quantized] = quantize(matrix,levels)
X = matrix(:); % Unroll into vector

% Round off to index
[~,X] = quantiz(X,levels(1:end-1),levels);

% quantized = ceil(matrix.*levels)./levels;

quantized = reshape(X,size(matrix)); % Reshape original matrix size
end