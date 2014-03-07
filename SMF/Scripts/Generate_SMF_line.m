% Output:
% SMF spatial matched filter
%
% Input:
% x_coord = x coordinate for line to generate 
% useCaseParams
% transducerType

function SMF = Generate_SMF_line(x_coord,resolution,useCaseParams,transducerType)

span = [0.01 0.05];
z = linspace(span(1),span(2), resolution(1))';
pos = num2cell([x_coord*ones(size(z)) zeros(size(z)) z],2);
SMF = cellfun( @(X) ...
    Generate_SMF_point(X,useCaseParams,transducerType),...
    pos,...
    'UniformOutput',false)

% SMF = cell(resolution(1),1);
% parfor i = 1:resolution(1)
%    SMF(i) = { Generate_SMF_point(pos{i},useCaseParams,transducerType) };
% end

end

function SMF = Generate_SMF_point(pos,useCaseParams,transducerType)

media.phantom_positions = pos;
media.phantom_amplitudes = 1;

f0 = 3e6;
fs = 120e6;
xmt_impulse_response = sin(2*pi*f0*(0:1/fs:2/f0))';
xmt_impulse_response = xmt_impulse_response.*hanning(max(size(xmt_impulse_response)));
rcv_impulse_response = xmt_impulse_response;

excitation = (sin(2*pi*f0*(0:1/fs:2/f0)))';
excitation = excitation.*hanning(max(size(excitation)));

SMF = Data_Acquisition('usecaseparams',useCaseParams, ...
                      'transducertype',transducerType, ...
                      'xmt_impulse_response', xmt_impulse_response, ...
                      'xmt_impulse_response_fs',fs, ...
                      'rcv_impulse_response', rcv_impulse_response, ...
                      'rcv_impulse_response_fs',fs, ...
                      'excitation_waveform', excitation, ... 
                      'excitation_fs',fs, ...
                      'symmetric','symmetric',...
                      'media',media);

end