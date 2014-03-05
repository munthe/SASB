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
%addpath(['./Scripts/bft3-beta-1-19/src'])
addpath(['..'])

usecase_filename = 'WirePhantom75MHzSlidingFilterOffNormalPulse';

savepath = loadpath;
savepath_figures = '../Figures8804/';
figure_format = '-depsc2';

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

%% Create Spat  ial Matched Filter

pos = [0 0 0.0150];
Generate_SMF_point(pos,useCaseParams,transducerType)



%% Define view parameters

useCaseParams.scanparams(1).windowtissueq.x_tismin = -0.04;
useCaseParams.scanparams(1).windowtissueq.y_tismin = 0;                
useCaseParams.scanparams(1).windowtissueq.x_tismax = 0.04;
useCaseParams.scanparams(1).windowtissueq.y_tismax = 0.10;
useCaseParams.scanparams(1).scantype = 2;
