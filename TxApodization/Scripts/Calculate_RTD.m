%%%%%%%%%%%%%%%%%%%%% Calculate Round Trip Delay %%%%%%%%%%%%%%%%%%%%%%%%%%
% Martin Hemmsen (ph.D)
% 27.06.08 (d.m.y)
% Ver. 01
%
% Use: This function is used to calculate the round trip delay, caused by
%      transducer specification and excitation pulse.
%
% Function: [FWHM, varargout] = Calculate_RTD(excitation,impulse_response)
%
% Input:
%   excitation : the excitation pulse
%   impulse_response : the impulse response
%   fs : sample frequency
%
% Output:
%   RT_delay: delay in seconds
%
% Needed matlab files:
%   * (mat) none
%   * (m) Resample.m
%
% Example of use:
%   fs=50e6;            % Sampling frequency [Hz]
%   xdc.f0 = 7e6;      	%  Center frequency           [Hz]
%   impulse_response = sin(2*pi*xdc.f0*(0:1/fs:2/xdc.f0));
%   impulse_response = impulse_response.*hanning(max(size(impulse_response)))';
%   excitation = sin(2*pi*xdc.f0*(0:1/fs:2/xdc.f0));
%   [RT_delay] = Calculate_RTD(excitation,impulse_response)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
function [RT_delay] = Calculate_RTD(excitation,xmt_impulse_response,rcv_impulse_response,fs)
%% Upsample to get a more precise estimate of the maximum
    fs_new = fs*10;
    excitation_new = resample(excitation,fs_new,fs);
    xmt_impulse_response_new = resample(xmt_impulse_response(:),fs_new,fs);
    rcv_impulse_response_new = resample(rcv_impulse_response(:),fs_new,fs);
    
%% Calculate Round trip delay due to transducer specification  
    % calculate delay at transmitter
    RTfunction_new = conv(excitation_new,xmt_impulse_response_new);
    % calculate delay at receiver
    RTfunction_new = conv(RTfunction_new,rcv_impulse_response_new);
    % calculate envelope of the round trip function
    RTfunction_new = abs(hilbert(RTfunction_new));
    % find delay at maximum response
    [J_new,RTfunction_new] = max(RTfunction_new);
    % adjust for matlabs array indexing
    RT_delay = (RTfunction_new-1)/fs_new;
    

