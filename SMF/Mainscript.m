%% INIT
% clc
% clear all
% close all

base_path = './';
loadpath = './';
base_path_usecase = './';

addpath(['./Scripts'])
% addpath(['./Scripts/Field'])
% addpath(['./Scripts/ScanConvert'])
%addpath(['./Scripts/bft3-beta-1-19/src'])
addpath(genpath('../lib'))

usecase_filename = 'WirePhantom75MHzSlidingFilterOffNormalPulse';

savepath = loadpath;
savepath_figures = '../Figures8804/';
figure_format = '-depsc2';

%% Initilize scanning settings
% init_parameters

%% CFUtools for measuring
CFUtools_init % init CFUtools
rect = [-40*ones(1,9); 10:10:90; 80*ones(1,9); 10*ones(1,9)]';


%% Read useCase file
useCaseParams = Read_Usecase(fullfile(base_path_usecase,[usecase_filename '.dat']));

%% Redefine transmit parameter

useCaseParams.scanparams(1).SubAppertureSize = 1;
useCaseParams.scanparams(1).scantype = 0; % 1d array imaging

transducerType = '8804';

% Matching layer not important for simulation, set to zero
useCaseParams.acmodparams(1).layerthickness1 = 0;
useCaseParams.acmodparams(1).layerthickness2 = 0;
useCaseParams.acmodparams(1).layerthickness3 = 0;

useCaseParams.bfxmitparams(1).xmitfnum = 2;
% 0 boxcar, 1 hamming, 2 gauss, 3 hanning, 4 blackman, 5 bartlett
useCaseParams.bfxmitparams(1).xmitapodishape = 0;
% useCaseParams.bfxmitparams(1).xmitapodigausswidth = 0.7; 
useCaseParams.bfxmitparams(1).xmitfocus = 0.02;
useCaseParams.bfxmitparams(1).xmitapodilevels = [];

%% Redefine receive parameter
useCaseParams.bfrcvparams(1).rcvapodilevels = [];

%% Create Spatial Matched Filter
resolution = [5 3];
tic
SMF = Generate_SMF(resolution,useCaseParams,transducerType);
toc

%% Generate scatter field
sca_x = ([ 0  0  0  0  0  0  0  0  0]./1000)';
sca_y = ([zeros(size(sca_x,1),1)]./1000);
sca_z = ([15 25 35 45 55 65 75 85 95]./1000)';
media.phantom_positions = [sca_x sca_y sca_z];
media.phantom_amplitudes = ones(size(media.phantom_positions,1),1);

%% Create first stage

f0 = 3e6;
fs = 120e6;
xmt_impulse_response = sin(2*pi*f0*(0:1/fs:2/f0))';
xmt_impulse_response = xmt_impulse_response.*hanning(max(size(xmt_impulse_response)));
rcv_impulse_response = xmt_impulse_response;

excitation = (sin(2*pi*f0*(0:1/fs:2/f0)))';
excitation = excitation.*hanning(max(size(excitation)));

RFdata = Data_Acquisition('usecaseparams',useCaseParams, ...
                      'transducertype',transducerType, ...
                      'xmt_impulse_response', xmt_impulse_response, ...
                      'xmt_impulse_response_fs',fs, ...
                      'rcv_impulse_response', rcv_impulse_response, ...
                      'rcv_impulse_response_fs',fs, ...
                      'excitation_waveform', excitation, ... 
                      'excitation_fs',fs, ...
                      'symmetric','symmetric',...
                      'media',media);
                  
                  
%% Create second stage image

image = Second_Stage_SMF(RFdata(1:3500,1:192),SMF, useCaseParams);


%% Plot filtered second stage image

setupDesc = [' Focus', num2str(useCaseParams.bfxmitparams(1).xmitfocus), ...
    ' F#',num2str(useCaseParams.bfxmitparams(1).xmitfnum) ...
    ];

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
plot_filtered_sasb_2 = Plot_Beamformed_Data('RFdata', image, ...
        'figure nr', 1, ...
        'data type','simulation',...
        'type of beamformation', 'sasb',...
        'stage','second',...
        'show tgc','dont show tgc',...
        'loadpath',savepath,...
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
% psf = cfu_get_psf_metrics('script',1,'fh',figure(fig_nr),'rect',rect(1,:));
% text(rect(1,1),rect(1,2),...
%     {['FWHM ' num2str(psf.fwhm_x) ' ' num2str(psf.fwhm_y)];...
%         ['Radius 20dB ' num2str(psf.radius20dB)]},...
%     'Color','red','FontSize',14,'VerticalAlignment','top');
% drawnow

print(figure(fig_nr),[savepath_figures 'SecondStageFiltered ' setupDesc],figure_format)


%% Plot
%imagesc(SMF(3,2).filter)
%colormap(gray)