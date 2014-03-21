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
% resolution = [5 3];
% tic
% SMF = Generate_SMF(resolution,useCaseParams,transducerType);
% toc

%% Generate scatter field
sca_z = ((15:10:95)/1000)';
sca_x = (zeros(size(sca_z))/1000);
sca_y = (zeros(size(sca_x))/1000);
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

SMFpath = '/usr/local/s103303/SMF_3000x43_crop60dB/';
resolution = [3000,43];
image = Second_Stage_SMF(RFdata,SMFpath,resolution,useCaseParams);
save(['./SMFimage' num2str(resolution(1)) 'x' num2str(resolution(1))], 'image','transducerType','useCaseParams','RFdata','resolution','SMFpath','media');

%% Plot filtered second stage image

% Calculate sampling frequency for secondstage image plotting
useCaseParams.bfrcvparams(1).smpfreq = size(image,1)*useCaseParams.scanparams(1).c_sound/(2*useCaseParams.scanparams(1).windowtissueq.y_tismax);

setupDesc = [' Focus', num2str(useCaseParams.bfxmitparams(1).xmitfocus), ...
    ' F#',num2str(useCaseParams.bfxmitparams(1).xmitfnum) ...
    ];

fig_nr = 2;
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
        'usecaseparams',useCaseParams,...
        'figure nr', fig_nr, ...
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

print(figure(fig_nr),[savepath_figures 'Filtered' setupDesc],figure_format)


%% Calculate point spread function

% % psf(1,9) = struct
% psf(setting,1) = cfu_get_psf_metrics('script',1,'fh',figure(fig_nr),'rect',rect(1,:));
% % psf(setting,1).setupDesc = setupDesc;
for i = 1:size(rect,1)
    psf(i) = cfu_get_psf_metrics('script',1,'fh',figure(fig_nr),'rect',rect(i,:));
%     psf(setting,i).setupDesc = setupDesc;
end


%% Plot psf as function of depth

fig_nr = 3;
figure(fig_nr);
depth = cell2mat({psf(:).y_coord})*1000;
plot( ...
    depth,cell2mat({psf(:).fwhm_x}),'o', ...
    depth,cell2mat({psf(:).fwhm_y}),'o', ...
    depth,cell2mat({psf(:).radius20dB}),'o' ...
    );
xlabel('Depth of point scatter [mm]');
ylabel('[mm]');
legend('FWHM_x','FWHM_y','Radius 20dB','Location','NorthWest');
% title(psf(setting,1).setupDesc);
prettyfig
print(figure(fig_nr),[savepath_figures 'PSF_SMF_' setupDesc],figure_format)


%% Plot
imagesc(SMF(3,2).filter)
colormap(gray)