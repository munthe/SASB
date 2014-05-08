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
savepath_figures = '../Figures_filtertest/';
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
% resolution = [2 2];
% tic
% SMF = Generate_SMF(resolution,useCaseParams,transducerType);
% toc

%% Generate scatter field
sca_z = ((15:10:95)/1000)';
sca_x = (zeros(size(sca_z))*15/1000);
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

%% Load filter
coord = [-19 75]/1000;
resolution = [5405,192];
line = round( (0.02+coord(1))/0.04 * resolution(2));
SMFpath = '/data/cfudata3/mah/Spatial_matched_filter/SMF_5405x192_f03MHz_c20dB/';
resolution = [5405,192];
load([SMFpath 'SMF_line_' num2str(line)], 'SMFline');

% Length and width of image [mm]
l =(useCaseParams.scanparams(1).windowtissueq.y_tismax-useCaseParams.scanparams(1).windowtissueq.y_tismin);
w =(useCaseParams.scanparams(1).windowtissueq.x_tismax-useCaseParams.scanparams(1).windowtissueq.x_tismin);

% Pixel depth of scatter
z = round( (coord(2)-useCaseParams.scanparams(1).windowtissueq.y_tismin)/l * resolution(1) );

% Define axis
ay = [(SMFline(z).index(1,1)/resolution(1))*l*1000 (SMFline(z).index(2,1)/resolution(1))*l*1000];
ax = [(SMFline(z).index(1,2)-resolution(2)/2)/resolution(2)*w*1000 (SMFline(z).index(2,2)-resolution(2)/2)/resolution(2)*w*1000];

% Choose filter
filter = zeros(size(RFdata));
filter(...
    SMFline(z).index(1,1):SMFline(z).index(2,1) ,...
    SMFline(z).index(1,2):SMFline(z).index(2,2) )...
    = SMFline(z).filter;

setupDesc = ['line_', num2str(line), 'Scatdepth', num2str(coord(2)*1000), ...
    ];

% Show filter
figure(3)
fig_nr = 3;
imagesc([useCaseParams.scanparams(1).windowtissueq.x_tismin*1000 useCaseParams.scanparams(1).windowtissueq.x_tismax*1000],...
    [useCaseParams.scanparams(1).windowtissueq.y_tismin*1000 useCaseParams.scanparams(1).windowtissueq.y_tismax*1000],filter)
colormap(linspecer)
colorbar
c = max(abs([min(filter(:)),max(filter(:))]));
caxis([-c c]);
xlabel('mm');
ylabel('mm');
prettyfig(8)
% Save
print(figure(fig_nr),[savepath_figures 'Filter_' setupDesc '.eps'],figure_format)

% Show cropped RFdata
RFdata_crop = RFdata(SMFline(z).index(1,1):SMFline(z).index(2,1) ,...
    SMFline(z).index(1,2):SMFline(z).index(2,2) );
figure(4)
fig_nr = 4;
imagesc([ax(1) ax(2)], [ay(1) ay(2)],RFdata_crop)
colormap(linspecer)
colorbar
c = max(abs([min(RFdata_crop(:)),max(RFdata_crop(:))]));
caxis([-c c]);
prettyfig(8)
xlabel('mm');
ylabel('mm');
% Save
print(figure(fig_nr),[savepath_figures 'RFdataCrop_' setupDesc '.eps'],figure_format)

% Show cropped filter
figure(5)
fig_nr = 5;
imagesc([ax(1) ax(2)], [ay(1) ay(2)],SMFline(z).filter)
colormap(linspecer)
colorbar
c = max(abs([min(SMFline(z).filter(:)),max(SMFline(z).filter(:))]));
caxis([-c c]);
prettyfig(8)
xlabel('mm');
ylabel('mm');
% Save
print(figure(fig_nr),[savepath_figures 'FilterCrop_' setupDesc '.eps'],figure_format)

% Show filtered data 
filtered = RFdata_crop.*SMFline(z).filter;
figure(6)
fig_nr = 6;
imagesc([ax(1) ax(2)], [ay(1) ay(2)],filtered)
colormap(linspecer)
colorbar
c = max(abs([min(filtered(:)),max(filtered(:))]));
caxis([-c c]);
prettyfig(8)
xlabel('mm');
ylabel('mm');
% Save
print(figure(fig_nr),[savepath_figures 'Filtered_' setupDesc '.eps'],figure_format)

%[m,i]=max(max(filter,[],1));
%fprintf('Max value (%e) at line %i, should be %i\n',m,i,(line-1)*size(RFdata,2)/resolution(2)+1)

%%
% a = zeros(192,1);
% b = zeros(192,1);
% for line = 1:192
%     SMFline=SMF(:,line);
% filter = zeros(size(RFdata));
% filter(...) ,...
%     SMFline(z).index(1,2):SMFline(z).index(2,2) )...
%     = SMFline(z).filter;
% [~,a(line)]=max(max(filter,[],1));
% b(line) = (line-1)*size(RFdata,2)/resolution(2)+1;
% end
% plot(a-b)
% %%
% fig_nr = 1;
% type = 'psf';
% switch(type)
%     case{'psf'}
%         compression = 'Linear';
%         dynamic_range = 60;
%     case{'B-mode'}
%         compression = 'BK_muLaw_RTSC_VP_DRC_SEL_0';
%         dynamic_range = 100;
% end
% clf(figure(fig_nr))
% plot_obj_sasb_2 = Plot_Beamformed_Data('RFdata', filter, ...
%         'usecaseparams',useCaseParams,...
%         'figure nr', fig_nr, ...
%         'data type','simulation',...
%         'type of beamformation', 'sasb',...
%         'stage','first',...
%         'show tgc','dont show tgc',...
%         'loadpath',[savepath ],...
%         'gain',[],...
%         'windowtissueq',useCaseParams.scanparams(1).windowtissueq,...
%         'compression',compression,...
%         'dynamic range',dynamic_range);
% axis image
% drawnow
% caxis([-60 0])
% set(gcf,'position',[  939   100   735   885])
% set(gca,'position',[ 80    70   640   800])
% drawnow

