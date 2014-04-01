
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

%% Load first stage
load cystimage_RFdata

%% Create second stage image

SMFpath = '/data/cfudata3/mah/Spatial_matched_filter/SMF_3000x192_crop60dB/';
resolution = [3000,192];
image = Second_Stage_SMF(RFdata,SMFpath,resolution,useCaseParams);
save(['./cystimage' num2str(resolution(1)) 'x' num2str(resolution(2))], 'image','transducerType','useCaseParams','RFdata','resolution','SMFpath','media');
