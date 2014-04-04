%% INIT
% clc
% clear all
% close all

base_path = './';
loadpath = '/data/cfudata3/mah/Spatial_matched_filter/cystphantom/';
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

%% Load RFdata

scanline = 1:192;
RFdata_collect = zeros(10000,length(scanline));
samples = zeros(length(scanline),1);
for i = 1:length(scanline);
    load([loadpath 'tmp/' 'cystRFdata_line' num2str(scanline(i))]);
    RFdata_collect(1:size(RFdata,1),i) = RFdata(:,scanline(i));
    samples(i) = size(RFdata,1);
end
RFdata = RFdata_collect(1:max(samples),:);

%% Save RFdata

str = [savepath 'cystRFdata' 'F' num2str(useCaseParams.bfxmitparams(1).xmitfnum) '_focus' num2str(useCaseParams.bfxmitparams(1).xmitfocus*1000) 'mm' 'f0' num2str(useCaseParams.bfxmitparams(1).xmitfreq/1e6) 'MHz'];
fprintf(['Saving RFdata to ' str '.mat ... ']);
save(str, 'transducerType','useCaseParams','RFdata','media');
fprintf('Saved.\n');