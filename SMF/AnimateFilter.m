%% INIT
% clc
% clear all
% close all

% Script to visualize filters

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
coord = [0 75]/1000;
resolution = [5405,192];
SMFpath = './';
load([SMFpath 'SMF_depth_' num2str(coord(2)*1000)], 'SMFdepth','useCaseParams');

AX = 1000*linspace(useCaseParams.scanparams(1).windowtissueq.x_tismin,useCaseParams.scanparams(1).windowtissueq.x_tismax,resolution(2));

for i=  resolution(2)/2:-1:1
    scanline = resolution(2) - i+1;
    SMFdepth(scanline).filter = fliplr(SMFdepth(i).filter);
    SMFdepth(scanline).index(:,1) = SMFdepth(i).index(:,1)
    SMFdepth(scanline).index(:,2) = [resolution(2)-SMFdepth(i).index(2,2)+1;resolution(2)-SMFdepth(i).index(1,2)+1];
end

writerObj = VideoWriter('filter.avi');
open(writerObj);

clear profile
for scanline = 1:resolution(2);

% Length and width of image [mm]
l =(useCaseParams.scanparams(1).windowtissueq.y_tismax-useCaseParams.scanparams(1).windowtissueq.y_tismin);
w =(useCaseParams.scanparams(1).windowtissueq.x_tismax-useCaseParams.scanparams(1).windowtissueq.x_tismin);

% Pixel depth of scatter
z = round( (coord(2)-useCaseParams.scanparams(1).windowtissueq.y_tismin)/l * resolution(1) );

% Define axis
ay = [(SMFdepth(scanline).index(1,1)/resolution(1))*l*1000 (SMFdepth(scanline).index(2,1)/resolution(1))*l*1000];
ax = [(SMFdepth(scanline).index(1,2)-resolution(2)/2)/resolution(2)*w*1000 (SMFdepth(scanline).index(2,2)-resolution(2)/2)/resolution(2)*w*1000];

% Choose filter
filter = zeros(size(RFdata));
filter(...
    SMFdepth(scanline).index(1,1):SMFdepth(scanline).index(2,1) ,...
    SMFdepth(scanline).index(1,2):SMFdepth(scanline).index(2,2) )...
    = SMFdepth(scanline).filter;

setupDesc = ['line_', num2str(line), 'Scatdepth', num2str(coord(2)*1000), ...
    ];

figure(1)

% Show cropped RFdata
subplot(231)
RFdata_crop = RFdata(SMFdepth(scanline).index(1,1):SMFdepth(scanline).index(2,1) ,...
    SMFdepth(scanline).index(1,2):SMFdepth(scanline).index(2,2) );
imagesc([ax(1) ax(2)], [ay(1) ay(2)],RFdata_crop)
colormap(linspecer)
% colorbar
c = max(abs([min(RFdata_crop(:)),max(RFdata_crop(:))]));
caxis([-c c]);
xlabel('mm');
ylabel('mm');
prettyfig

% Show cropped filter
subplot(232)
imagesc([ax(1) ax(2)], [ay(1) ay(2)],SMFdepth(scanline).filter)
colormap(linspecer)
% colorbar
c = max(abs([min(SMFdepth(scanline).filter(:)),max(SMFdepth(scanline).filter(:))]));
caxis([-c c]);
xlabel('mm');
ylabel('mm');
prettyfig

% Show filtered data 
subplot(233)
filtered = RFdata_crop.*SMFdepth(scanline).filter;
imagesc([ax(1) ax(2)], [ay(1) ay(2)],filtered)
colormap(linspecer)
% colorbar
c = max(abs([min(filtered(:)),max(filtered(:))]));
caxis([-c c]);
xlabel('mm');
ylabel('mm');
prettyfig

subplot(212)
profile(scanline) = sum(filtered(:));
scatter(AX(scanline),profile(scanline))
plot(AX(1:scanline),profile,'LineWidth',2)
axis([AX(1) AX(end) -2.1116e-45 6.1314e-44])

prettyfig
frame = getframe(gcf);
writeVideo(writerObj,frame);

end

close(writerObj);

