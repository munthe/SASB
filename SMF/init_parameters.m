% Read useCase file and initilize parameters for simulations.

base_path_usecase = './';
usecase_filename = 'WirePhantom75MHzSlidingFilterOffNormalPulse';
transducerType = '8804';


%% Read useCase file
useCaseParams = Read_Usecase(fullfile(base_path_usecase,[usecase_filename '.dat']));

%% Redefine transmit parameters
useCaseParams.scanparams(1).SubAppertureSize = 1;
useCaseParams.scanparams(1).scantype = 0; % 1d array imaging

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

%% Redefine receive parameters
useCaseParams.bfrcvparams(1).rcvapodilevels = [];

%% Define view parameters
useCaseParams.scanparams(1).windowtissueq.x_tismin = -0.02;
useCaseParams.scanparams(1).windowtissueq.y_tismin = 0;                
useCaseParams.scanparams(1).windowtissueq.x_tismax = 0.02;
useCaseParams.scanparams(1).windowtissueq.y_tismax = 0.10;
useCaseParams.scanparams(1).scantype = 2;
