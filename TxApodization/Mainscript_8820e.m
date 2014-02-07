%% INIT
% clc
% clear all
% close all

base_path = './';
loadpath = './';
base_path_usecase = './';

addpath(['./Scripts'])
addpath(['./Scripts/Field'])
addpath(['./Scripts/ScanConvert'])
addpath(['./Scripts/bft3-beta-1-19/src'])

usecase_filename = 'Pre_SASB_liver_2_1_p';

savepath = loadpath;
savepath_figures = './';

%% Read useCase file
useCaseParams = Read_Usecase(fullfile(base_path_usecase,[usecase_filename '.dat']));
useCaseParams.bfrcvparams(1) = useCaseParams.bfrcvparams(11);

%% Generate scatter field

sca_x = ([ 0  0  0  0  0  0  0  0  0]./1000)';
sca_y = ([zeros(size(sca_x,1),1)]./1000);
sca_z = ([15 25 35 45 55 65 75 85 95]./1000)';
% sca_x = 0;
% sca_y = 0;
% sca_z = 0.05;
media.phantom_positions = [sca_x sca_y sca_z];
media.phantom_amplitudes = ones(size(media.phantom_positions,1),1);

%% Redefine transmit parameter

useCaseParams.scanparams(1).SubAppertureSize = 1;
useCaseParams.scanparams(1).scantype = 0; % 1d array imaging

transducerType = '8820e';

% Matching layer not important for simulation, set to zero
useCaseParams.acmodparams(1).layerthickness1 = 0;
useCaseParams.acmodparams(1).layerthickness2 = 0;
useCaseParams.acmodparams(1).layerthickness3 = 0;

useCaseParams.bfxmitparams(1).xmitfnum = 1;
% 0 boxcar, 1 hamming, 2 gauss, 3 hanning, 4 blackman, 5 bartlett
useCaseParams.bfxmitparams(1).xmitapodishape = 1;
useCaseParams.bfxmitparams(1).xmitapodigausswidth = 0.7; 
useCaseParams.bfxmitparams(1).xmitfocus = 0.02;
% quantization of apodization round of to levels in vector, if empty not
% rounding
useCaseParams.bfxmitparams(1).xmitapodilevels = [];

%% Set the impulse response and excitation

f0 = 3e6;
fs = 120e6;
xmt_impulse_response = sin(2*pi*f0*(0:1/fs:2/f0))';
xmt_impulse_response = xmt_impulse_response.*hanning(max(size(xmt_impulse_response)));
rcv_impulse_response = xmt_impulse_response;

excitation = (sin(2*pi*f0*(0:1/fs:2/f0)))';
excitation = excitation.*hanning(max(size(excitation)));

%% Generate Data - Beamform first stage

Data_AcquisitionVer3('usecaseparams',useCaseParams, ...
                      'transducertype',transducerType, ...
                      'xmt_impulse_response', xmt_impulse_response, ...
                      'xmt_impulse_response_fs',fs, ...
                      'rcv_impulse_response', rcv_impulse_response, ...
                      'rcv_impulse_response_fs',fs, ...
                      'excitation_waveform', excitation, ... 
                      'excitation_fs',fs, ...
                      'symmetric','symmetric',...
                      'savepath',savepath, ...
                      'media',media)

%% Define second stage beamforming parameters

useCaseParams.bfrcvparams(1).rcvfnum = useCaseParams.bfxmitparams(1).xmitfnum;
useCaseParams.bfrcvparams(1).rcvapodishape = 1;
useCaseParams.bfrcvparams(1).rcvapodigausswidth = 0.2;
useCaseParams.bfrcvparams(1).rcvfocus = useCaseParams.bfxmitparams(1).xmitfocus;

useCaseParams.scanparams(1).startdepthq = 0.0;
useCaseParams.scanparams(1).stopdepthq = 0.12;

useCaseParams.scanparams(1).startlinenumq = 0;
useCaseParams.scanparams(1).stoplinenumq = 269;

useCaseParams.scanparams(1).scanareadscr.startlineorigin.y = 0;
useCaseParams.scanparams(1).scanareadscr.startlineorigin.x = 0.04;
useCaseParams.scanparams(1).scanareadscr.startlineangle = pi/2;

useCaseParams.scanparams(1).scanareadscr.stoplineorigin.y = 0;
useCaseParams.scanparams(1).scanareadscr.stoplineorigin.x = -0.04;
useCaseParams.scanparams(1).scanareadscr.stoplineangle =  pi/2;

%useCaseParams.scanparams(1).steeringangle = 0;
%useCaseParams.scanparams(1).phasedangle = 30/180*pi;

%% Beamform 2nd stage

tic
BeamformationVer4('beamformation_method','sasbsim_new', ...
                  'usecase',useCaseParams, ...
                  'loadpath', savepath, ...
                  'symmetric','symmetric',...
                  'savepath',[savepath],...
                  'generate_scaling',false)
toc

%% Define view parameters

useCaseParams.scanparams(1).windowtissueq.x_tismin = -0.04;
useCaseParams.scanparams(1).windowtissueq.y_tismin = 0;                
useCaseParams.scanparams(1).windowtissueq.x_tismax = 0.04;
useCaseParams.scanparams(1).windowtissueq.y_tismax = 0.10;
useCaseParams.scanparams(1).scantype = 2;

%% Plot first stage image

fig_nr = 1;
type = 'psf';
switch(type)
    case{'psf'}
        compression = 'Linear';
        dynamic_range = 60;
    case{'B-mode'}
        compression = 'BK_muLaw_RTSC_VP_DRC_SEL_0';
        dynamic_range = 100;
end
clf(figure(fig_nr))
plot_obj_sasb_2 = Plot_Beamformed_Data_Ver2(...
        'figure nr', fig_nr, ...
        'data type','simulation',...
        'type of beamformation', 'sasb',...
        'stage','first',...
        'show tgc','dont show tgc',...
        'loadpath',[savepath ],...
        'gain',[],...
        'windowtissueq',useCaseParams.scanparams(1).windowtissueq,...
        'compression',compression,...
        'dynamic range',dynamic_range); 
axis image   
drawnow
caxis([-60 0])
set(gcf,'position',[  939   100   735   885])
set(gca,'position',[ 80    70   640   800])
drawnow

%% Plot second stage image

fig_nr = 8;
type = 'psf';
switch(type)
    case{'psf'}
        compression = 'Linear';
        dynamic_range = 60;
    case{'B-mode'}
        compression = 'BK_muLaw_RTSC_VP_DRC_SEL_0';
        dynamic_range = 100;
end
plot_obj_sasb_2 = Plot_Beamformed_Data_Ver2(...
        'figure nr', 2, ...
        'data type','simulation',...
        'type of beamformation', 'sasb',...
        'stage','second',...
        'show tgc','dont show tgc',...
        'loadpath',[savepath ],...
        'gain',[],...
        'windowtissueq',useCaseParams.scanparams(1).windowtissueq,...
        'compression',compression,...
        'dynamic range',dynamic_range); 
    
axis image   
drawnow
caxis([-60 0])
set(gcf,'position',[  939   100   735   885])
set(gca,'position',[ 80    70   640   800])
drawnow
