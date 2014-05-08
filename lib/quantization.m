% Function for quantization of the apodization vector
%
% Inputs:
% Matrix to be quantized
% Quantization levels in an ascending vector
%
% Outputs:
% Quantized matrix

function [quantized] = quantization(matrix,levels)
quantized = zeros(size(matrix));
for i = 1:length(levels)-1;
    % Mid value for rounding
    mid_value = (levels(i)+levels(i+1))/2;
    
    % Set values greater than the current level and up to the midvalue to
    % current level
    quantized ( levels(i)<matrix & matrix<mid_value ) = levels(i);
    % Set values greater than midvalue and up to next level og next level
    quantized ( mid_value<=matrix & matrix<levels(i+1) ) = levels(i+1);
end
end