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

function SMF = Generate_SMF(resolution,useCaseParams,transducerType)

span = [-10,10]/1000;
x_coord = linspace(span(1),span(2), resolution(2));

SMF = arrayfun( @(x) ...
    Generate_SMF_line(x,resolution,useCaseParams,transducerType),...
    x_coord,...
    'UniformOutput',false);
SMF = [SMF{:}];

% SMF = cell(resolution);
% parfor i = 1:resolution(2)
%     SMF(:,i) = Generate_SMF_line(x_coord(i),useCaseParams,transducerType);
% end

end