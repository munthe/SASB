% Output:
% SMF spatial matched filter
%
% Input:
% x_coord = x coordinate for line to generate 
% useCaseParams
% transducerType

function SMF = Generate_SMF_line(x_coord,useCaseParams,transducerType)

span = [0.01 0.05];
num = 5;
z = linspace(span(1),span(2), num)';
pos = num2cell([x_coord*ones(size(z)) zeros(size(z)) z],2);
SMF = cellfun( @(X) ...
    Generate_SMF_point(X,useCaseParams,transducerType),...
    pos,...
    'UniformOutput',false);

end