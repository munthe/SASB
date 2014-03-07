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

savepath = loadpath;
savepath_figures = '../Figures8804/';
figure_format = '-depsc2';

%% Initilize scanning settings
init_parameters

%% CFUtools for measuring
CFUtools_init % init CFUtools
rect = [-40*ones(1,9); 10:10:90; 80*ones(1,9); 10*ones(1,9)]';

%% Create Spatial Matched Filter
resolution = [5 3];
tic
SMF = Generate_SMF(resolution,useCaseParams,transducerType);
toc

%% Plot
imagesc(SMF{1,3})
colormap(gray)