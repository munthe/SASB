function [samples_mf cutoff] = cfu_mf (samples, waveform, shape)

%
% SAMPLES_MF [CUTOFF] = CFU_MF(samples, waveform [,shape])
%
% The input samples are converted to double, the DC-offset is removed from each channel, and then 
% the matched-filtering is carried out.
%
% Before filtering the waveform amplitude is normalized to have unit area.
%
%
% EXAMPLE:
%[samples_mf cut_idx] = mf(samples, xmt.wf, 'full');% match filtering
%rf_data = hilbert(samples_mf); % IQ beamforming
%rf_data=rf_data(1+cut_idx:end-cut_idx,:); %remove extra samples from convolution
%
% 2012-05-25, V1.0, Init version, MFR
% 2012-09-07, v1.1, Help text expanded. MFR
% 2012-09-22, v1.2, Shape of conv2 now an input argument.
% 2012-09-27, v1.3, cutoff index is now calculated when shape='full'.
% 2013-09-17, v1.4, Renamed to cfu_mf. MFR
%

if  nargin < 3,    shape = 'same'; end
if ~ischar(shape), error('''shape'' must be a string.'); end
    


% kernel
kernel  = fliplr(waveform)';
kernel  = kernel/sum(abs(kernel));
% The following forces 'samples' to be copied, not just referenced. For speed, maybe
% the next two lines should be left out..
samples = double(samples);                                     % parse to double
samples = samples - repmat(mean(samples),[size(samples,1) 1]); % remove DC from signal
% match filtering
samples_mf = conv2(samples, kernel, shape);


% set second output
if nargout > 1
    if strcmp(shape,'full')
        cutoff = floor(length(waveform)/2);
    else
        cutoff = 0;
    end
end


