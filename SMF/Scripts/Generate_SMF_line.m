% Output:
% SMF spatial matched filter
%
% Input:
% x_coord = x coordinate for line to generate 
% useCaseParams
% transducerType

function SMF = Generate_SMF_line(x_coord,useCaseParams,transducerType)

span = [0.015 0.055];
num = 5;

z_coord = linspace(span(1),span(2), num);
SMF = {};
for i = 1:num
    pos = [x_coord 0 z_coord(i)];
    SMF(i) = { Generate_SMF_point(pos,useCaseParams,transducerType) };
end

end