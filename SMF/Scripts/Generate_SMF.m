% Function to generate a spatial matched filter. For use with second stage
% spatial matched filter beamforming.
% 
% Output:
% SMF - spatial matched filter
%
% Input:
% 
% useCaseParams
% transducerType

function SMF = Generate_SMF(useCaseParams,transducerType)

span = [-10,10]/1000;
num = 3;
x_coord = linspace(span(1),span(2), num);
SMF = arrayfun( @(x) ...
    Generate_SMF_line(x,useCaseParams,transducerType),...
    x_coord,...
    'UniformOutput',false);
SMF = [SMF{:}];

end