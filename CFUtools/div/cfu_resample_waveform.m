function wf_out = cfu_resample_waveform (wf, fs_input, fs_output)
% 
% RESAMPLE_WAVEFORM
% This function resamples a given waveform to a new sampling frequency.
%
% wf_out = cfu_resample_waveform (wf, fs_input, fs_output)
%
% By MOFI, 2013
% 2013-02-22, Init Version
% 2013-09-17, Renamed to cfu_resample_waveform. 
%

if fs_input == 0, error('fs_input may not be zero.'); end
    
[P,Q] = rat(fs_output/fs_input);
wf_out = resample(wf,P,Q);

