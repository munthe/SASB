% Function to quantize the apodization vector
%
% Inputs:
%
% Outputs:
%

function [quantized] = quantize(matrix,levels)

quantized = ceil(matrix.*levels)./levels;



end