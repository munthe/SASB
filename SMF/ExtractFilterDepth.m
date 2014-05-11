%% INIT
% clc
% clear all
% close all

% Script to extract filter at defined depth from all filter lines

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
savepath_figures = '../Figures_filtertest/';
figure_format = '-depsc2';
setupfigdefault

%% Load filter
coord = [0 75]/1000;
resolution = [5405,192];
line = round( (0.02+coord(1))/0.04 * resolution(2));
SMFpath = './filter/SMF_5405x192_f03MHz_c20dB/';
resolution = [5405,192];
load([SMFpath 'SMF_line_1'], 'useCaseParams')
% Length and width of image [mm]
l =(useCaseParams.scanparams(1).windowtissueq.y_tismax-useCaseParams.scanparams(1).windowtissueq.y_tismin);
w =(useCaseParams.scanparams(1).windowtissueq.x_tismax-useCaseParams.scanparams(1).windowtissueq.x_tismin);
% Pixel depth of scatter
z = round( (coord(2)-useCaseParams.scanparams(1).windowtissueq.y_tismin)/l * resolution(1) );
% Load lines and save specified depth
for l = 1:resolution(2)/2
    load([SMFpath 'SMF_line_' num2str(l)], 'SMFline');
    SMFdepth(l) = SMFline(z,1);
    clear SMFline
end
save([savepath 'SMF_depth_' num2str(coord(2)*1000)],'SMFdepth','useCaseParams');
