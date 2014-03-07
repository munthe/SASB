% Generate one line of the spatial matched filter, for use with SMF_slurm.m

addpath('../cluster');
addpath(genpath('../lib'));
addpath('./Scripts');

% init_parameters;

%% CFUtools for measuring
CFUtools_init % init CFUtools
rect = [-40*ones(1,9); 10:10:90; 80*ones(1,9); 10*ones(1,9)]';

%% Read useCase file
base_path_usecase = './';
usecase_filename = 'WirePhantom75MHzSlidingFilterOffNormalPulse';
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


savepath = '/data/cfudata3/mah/Spatial_matched_filter/';
tmpdir = 'tmp_SMF/';

line = par.line;
resolution = [5,3];

x_coord = linspace(useCaseParams.scanparams(1).windowtissueq.x_tismax,...
                   useCaseParams.scanparams(1).windowtissueq.x_tismin,...
                   resolution(2));

SMFline = Generate_SMF_line(x_coord(line),resolution(1),useCaseParams,transducerType);

%% Saving SMF
save([savepath tmpdir 'SMF_line_' num2str(line)],'SMFline')
