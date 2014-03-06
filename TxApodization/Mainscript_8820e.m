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

usecase_filename = 'Pre_SASB_liver_2_1_p';

savepath = loadpath;
savepath_figures = '../Figures/';
figure_format = '-depsc2';
setupfigdefault

%% CFUtools for measuring
CFUtools_init % init CFUtools
rect = [-40*ones(1,9); 10:10:90; 80*ones(1,9); 10*ones(1,9)]';

%% Generate scatter field
sca_x = ([ 0  0  0  0  0  0  0  0  0]./1000)';
sca_y = ([zeros(size(sca_x,1),1)]./1000);
sca_z = ([15 25 35 45 55 65 75 85 95]./1000)';
media.phantom_positions = [sca_x sca_y sca_z];
media.phantom_amplitudes = ones(size(media.phantom_positions,1),1);

%% Loop through different settings
exp = [0,2.^[0:14]];
apodizations = { ...
    [0,1],              1,  'Apo boxcar'    ; ...
    [],                 2,  'Apo inf'       ; ...  
    linspace(0,1,16),   3,  'Apo lin 16'    ; ...
    linspace(0,1,8),    3,  'Apo lin 8'     ; ...
    linspace(0,1,6),    3,  'Apo lin 6'     ; ...
    linspace(0,1,5),    3,  'Apo lin 5'     ; ...
    linspace(0,1,4),    3,  'Apo lin 4'     ; ...
    exp(1:6)/exp(6),    4,  'Apo exp 6'     ; ...
    exp(1:5)/exp(5),    4,  'Apo exp 5'     ; ...
    exp(1:4)/exp(4),    4,  'Apo exp 4'     ; ...
    [0,0.9,1],          5,  'Apo 0,0.9,1'   ; ...
    [0,0.8,1],          5,  'Apo 0,0.8,1'   ; ...
    [0,0.7,1],          5,  'Apo 0,0.7,1'   ; ...
    [0,0.6,1],          5,  'Apo 0,0.6,1'   ; ...
    [0,0.5,1],          5,  'Apo 0,0.5,1'   ; ...
    [0,0.4,1],          5,  'Apo 0,0.4,1'   ; ...
    [0,0.3,1],          5,  'Apo 0,0.3,1'   ; ...
    [0,0.2,1],          5,  'Apo 0,0.2,1'   ; ...
    [0,0.1,1],          5,  'Apo 0,0.1,1'   ; ...
    };
% apoDesc = {'Apo inf'};
% apodizations = {[]};
setting = 1;
%%
for setting = 1:size(apodizations,1) %rest of script
fprintf('Beginning simulation no %i\n',setting);

%% Read useCase file
useCaseParams = Read_Usecase(fullfile(base_path_usecase,[usecase_filename '.dat']));
useCaseParams.bfrcvparams(1) = useCaseParams.bfrcvparams(11);

%% Redefine transmit parameter

useCaseParams.scanparams(1).SubAppertureSize = 1;
useCaseParams.scanparams(1).scantype = 0; % 1d array imaging

transducerType = '8820e';

% Matching layer not important for simulation, set to zero
useCaseParams.acmodparams(1).layerthickness1 = 0;
useCaseParams.acmodparams(1).layerthickness2 = 0;
useCaseParams.acmodparams(1).layerthickness3 = 0;

useCaseParams.bfxmitparams(1).xmitfnum = 1.5;
% 0 boxcar, 1 hamming, 2 gauss, 3 hanning, 4 blackman, 5 bartlett
useCaseParams.bfxmitparams(1).xmitapodishape = 1;
% useCaseParams.bfxmitparams(1).xmitapodigausswidth = 0.7; 
useCaseParams.bfxmitparams(1).xmitfocus = 0.02;

switch( length(apodizations{setting,1}) )
    case(2) % use boxcar insted of rounding off
        useCaseParams.bfxmitparams(1).xmitapodishape = 0;
        useCaseParams.bfxmitparams(1).xmitapodilevels = [];
    otherwise % round of to levels in vector
        useCaseParams.bfxmitparams(1).xmitapodilevels = apodizations{setting,1};
end

%% Redefine receive parameter
useCaseParams.bfrcvparams(1).rcvapodilevels = [];

%% Set the impulse response and excitation

f0 = 3e6;
fs = 120e6;
xmt_impulse_response = sin(2*pi*f0*(0:1/fs:2/f0))';
xmt_impulse_response = xmt_impulse_response.*hanning(max(size(xmt_impulse_response)));
rcv_impulse_response = xmt_impulse_response;

excitation = (sin(2*pi*f0*(0:1/fs:2/f0)))';
excitation = excitation.*hanning(max(size(excitation)));

%% Generate Data - Beamform first stage

Data_AcquisitionVer3('usecaseparams',useCaseParams, ...
                      'transducertype',transducerType, ...
                      'xmt_impulse_response', xmt_impulse_response, ...
                      'xmt_impulse_response_fs',fs, ...
                      'rcv_impulse_response', rcv_impulse_response, ...
                      'rcv_impulse_response_fs',fs, ...
                      'excitation_waveform', excitation, ... 
                      'excitation_fs',fs, ...
                      'symmetric','symmetric',...
                      'savepath',savepath, ...
                      'media',media)

%% Define second stage beamforming parameters

useCaseParams.bfrcvparams(1).rcvfnum = useCaseParams.bfxmitparams(1).xmitfnum;
useCaseParams.bfrcvparams(1).rcvapodishape = 1;
useCaseParams.bfrcvparams(1).rcvapodigausswidth = 0.2;
useCaseParams.bfrcvparams(1).rcvfocus = useCaseParams.bfxmitparams(1).xmitfocus;

useCaseParams.scanparams(1).startdepthq = 0.0;
useCaseParams.scanparams(1).stopdepthq = 0.12;

useCaseParams.scanparams(1).startlinenumq = 0;
useCaseParams.scanparams(1).stoplinenumq = 268;

useCaseParams.scanparams(1).scanareadscr.startlineorigin.y = 0;
useCaseParams.scanparams(1).scanareadscr.startlineorigin.x = 0.04;
useCaseParams.scanparams(1).scanareadscr.startlineangle = pi/2;

useCaseParams.scanparams(1).scanareadscr.stoplineorigin.y = 0;
useCaseParams.scanparams(1).scanareadscr.stoplineorigin.x = -0.04;
useCaseParams.scanparams(1).scanareadscr.stoplineangle =  pi/2;

%useCaseParams.scanparams(1).steeringangle = 0;
%useCaseParams.scanparams(1).phasedangle = 30/180*pi;

%% Define view parameters

useCaseParams.scanparams(1).windowtissueq.x_tismin = -0.04;
useCaseParams.scanparams(1).windowtissueq.y_tismin = 0;                
useCaseParams.scanparams(1).windowtissueq.x_tismax = 0.04;
useCaseParams.scanparams(1).windowtissueq.y_tismax = 0.10;
useCaseParams.scanparams(1).scantype = 2;

%% Beamform 2nd stage

tic
BeamformationVer4('beamformation_method','sasbsim_new', ...
                  'usecase',useCaseParams, ...
                  'loadpath', savepath, ...
                  'symmetric','symmetric',...
                  'savepath',[savepath],...
                  'generate_scaling',false)
toc

%% Plot first stage image
setupDesc = [apodizations{setting,3}, ...
     ' Focus', num2str(useCaseParams.bfxmitparams(1).xmitfocus), ...
    ' F#',num2str(useCaseParams.bfxmitparams(1).xmitfnum) ...
    ];

fig_nr = 1;
type = 'psf';
switch(type)
    case{'psf'}
        compression = 'Linear';
        dynamic_range = 60;
    case{'B-mode'}
        compression = 'BK_muLaw_RTSC_VP_DRC_SEL_0';
        dynamic_range = 100;
end
clf(figure(fig_nr))
plot_obj_sasb_2 = Plot_Beamformed_Data_Ver2(...
        'figure nr', fig_nr, ...
        'data type','simulation',...
        'type of beamformation', 'sasb',...
        'stage','first',...
        'show tgc','dont show tgc',...
        'loadpath',[savepath ],...
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
prettyfig
print(figure(fig_nr),[savepath_figures 'FirstStage ' setupDesc],figure_format)


%% Plot second stage image

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
plot_obj_sasb_2 = Plot_Beamformed_Data_Ver2(...
        'figure nr', 2, ...
        'data type','simulation',...
        'type of beamformation', 'sasb',...
        'stage','second',...
        'show tgc','dont show tgc',...
        'loadpath',[savepath ],...
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
prettyfig
print(figure(fig_nr),[savepath_figures 'SecondStage ' setupDesc],figure_format)


%% Calculate point spread function

% % psf(1,9) = struct
% psf(setting,1) = cfu_get_psf_metrics('script',1,'fh',figure(fig_nr),'rect',rect(1,:));
% % psf(setting,1).setupDesc = setupDesc;
for i = 1:size(rect,1)
    psf(setting,i) = cfu_get_psf_metrics('script',1,'fh',figure(fig_nr),'rect',rect(i,:));
%     psf(setting,i).setupDesc = setupDesc;
end


%% Plot psf as function of depth

fig_nr = 3;
figure(fig_nr);
depth = cell2mat({psf(setting,:).y_coord})*1000;
plot( ...
    depth,cell2mat({psf(setting,:).fwhm_x}),'o', ...
    depth,cell2mat({psf(setting,:).fwhm_y}),'o', ...
    depth,cell2mat({psf(setting,:).radius20dB}),'o' ...
    );
xlabel('Depth of point scatter [mm]');
ylabel('[mm]');
legend('FWHM_x','FWHM_y','Radius 20dB','Location','NorthWest');
% title(psf(setting,1).setupDesc);
prettyfig
print(figure(fig_nr),[savepath_figures 'PSF SecondStage ' setupDesc],figure_format)

end


%% Plot comparrison of psf between simulations

setupDesc = [' Focus', num2str(useCaseParams.bfxmitparams(1).xmitfocus), ...
    ' F#',num2str(useCaseParams.bfxmitparams(1).xmitfnum)];
groups = cell2mat(apodizations(:,2));
c = linspecer(max(groups));
boxplotcolors = c(groups,:);

fig_nr=4;
figure(fig_nr);
fwhm_x = reshape(cell2mat({psf(:).fwhm_x}),size(psf));
fwhm_x_norm = fwhm_x(2:end,:)./repmat(fwhm_x(1,:),size(psf,1)-1,1);
boxplot(fwhm_x_norm',apodizations(2:end,3),'plotstyle','compact','colors',boxplotcolors(2:end,:));
prettyfig
print(figure(fig_nr),[savepath_figures 'Comparrison of PSF fwhm_x ' setupDesc],figure_format)

fig_nr=5;
figure(fig_nr);
radius20dB = reshape(cell2mat({psf(:).radius12dB}),size(psf));
radius20dB_norm = radius20dB(2:end,:)./repmat(radius20dB(1,:),size(psf,1)-1,1);
boxplot(radius20dB_norm',apodizations(2:end,3),'plotstyle','compact','colors',boxplotcolors(2:end,:));
prettyfig
print(figure(fig_nr),[savepath_figures 'Comparrison of PSF radius 20dB ' setupDesc],figure_format)

%% Save

save(strcat('Transducer ',transducerType, setupDesc))