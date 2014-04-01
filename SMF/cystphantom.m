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
% addpath(['./Scripts/bft3-beta-1-19/src'])
addpath(genpath('../lib'))

usecase_filename = 'WirePhantom75MHzSlidingFilterOffNormalPulse';

savepath = loadpath;
savepath_figures = '../Figures_SMF/';
figure_format = '-depsc2';
setupfigdefault

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

%% Define view parameters
useCaseParams.scanparams(1).windowtissueq.x_tismin = -0.02;
useCaseParams.scanparams(1).windowtissueq.y_tismin = 0.001;                
useCaseParams.scanparams(1).windowtissueq.x_tismax = 0.02;
useCaseParams.scanparams(1).windowtissueq.y_tismax = 0.101;
% useCaseParams.scanparams(1).scantype = 2;


%% Create Spatial Matched Filter
% resolution = [1000 192];
% tic
% SMF = Generate_SMF(resolution,useCaseParams,transducerType);
% toc

%% Generate scatter field
medium.x    = [-0.02 0.02]; % x-limit of medium in meters
medium.y    = [-0.005 0.005]; % y-limit of medium in meters
medium.z    = [0.04 0.09]; % z-limit of medium in meters
medium.dens = 0.5/(0.5*0.7); % average number of scatterer per mm^3

cyst.r      = [0.005 0.005 0.005]; % radius in the x, y and z-dimension
cyst.c      = [-0.01 0 0.075]; % center coordinate of the cyst 
cyst.amp    = 0; % amplitude of cyst (set to 0 for anechoic, 1 for same amplitude as medium)

[media.phantom_positions, media.phantom_amplitudes] = cfu_cyst_phantom(medium, cyst);

%% Create first stage

f0 = 3e6;
fs = 120e6;
xmt_impulse_response = sin(2*pi*f0*(0:1/fs:2/f0))';
xmt_impulse_response = xmt_impulse_response.*hanning(max(size(xmt_impulse_response)));
rcv_impulse_response = xmt_impulse_response;

excitation = (sin(2*pi*f0*(0:1/fs:2/f0)))';
excitation = excitation.*hanning(max(size(excitation)));

lines = [80,100];
RFdata_lines = arrayfun( ...
    @(l)(Data_Acquisition('usecaseparams',useCaseParams, ...
        'transducertype',transducerType, ...
        'xmt_impulse_response', xmt_impulse_response, ...
        'xmt_impulse_response_fs',fs, ...
        'rcv_impulse_response', rcv_impulse_response, ...
        'rcv_impulse_response_fs',fs, ...
        'excitation_waveform', excitation, ... 
        'excitation_fs',fs, ...
        'symmetric','symmetric',...
        'media',media,...
        'scanlines',l)),...
    lines,...
    'UniformOutput',false);
% Find max depth of scanlines:
depth = max(arrayfun(@(x)(size(RFdata{x},1)),1:size(RFdata,2)));
RFdata = zeros(depth,lines);

%% Save
str = ['./cystimage_RFdata'];
fprintf(['Saving RFdata to ' str '.mat ... ']);
save(str, 'transducerType','useCaseParams','RFdata','media');
fprintf('Saved.\n');
                  
%% Create second stage image

SMFpath = '/data/cfudata3/mah/Spatial_matched_filter/SMF_3000x192_crop60dB/';
resolution = [3000,192];
image = Second_Stage_SMF(RFdata,SMFpath,resolution,useCaseParams);

%% Save
str = ['./cystimage' num2str(resolution(1)) 'x' num2str(resolution(2))];
fprintf(['Saving RFdata and SMF filtered image to ' str '.mat ... ']);
save(str, 'image','transducerType','useCaseParams','RFdata','resolution','SMFpath','media');
fprintf('Saved.\n');
