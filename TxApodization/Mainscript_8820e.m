%% INIT
% clc
% clear all
% close all


base_path = './';
loadpath = './';
base_path_usecase = './';


%addpath(genpath(fullfile(base_path,'Matlab code','Toolbox','Display_functions')))
addpath(genpath(fullfile(base_path,'Scripts')))
addpath(['./Scripts/Field'])
addpath(genpath(fullfile(base_path,'Scripts','bft3-beta-1-24','src')))


usecase_filename = 'Pre_SASB_liver_2_1_p';

savepath = loadpath;
savepath_figures = './';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define Transmit parameters
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
useCaseParams = Read_Usecase(fullfile(base_path_usecase,[usecase_filename '.dat']));
useCaseParams.bfrcvparams(1) = useCaseParams.bfrcvparams(11);


useCaseParams.scanparams(1).SubAppertureSize = 1;
useCaseParams.scanparams(1).scantype = 0;

transducerType = '8820e';

useCaseParams.acmodparams(1).layerthickness1 = 0;
useCaseParams.acmodparams(1).layerthickness2 = 0;
useCaseParams.acmodparams(1).layerthickness3 = 0;
     

useCaseParams.bfxmitparams(1).xmitfnum = 2;
useCaseParams.bfxmitparams(1).xmitapodishape = 1;
useCaseParams.bfxmitparams(1).xmitapodigausswidth = 0.2;
useCaseParams.bfxmitparams(1).xmitfocus = 0.04;


%  Set the impulse response and excitation
f0 = 3e6;
fs = 120e6;
xmt_impulse_response = sin(2*pi*f0*(0:1/fs:2/f0))';
xmt_impulse_response = xmt_impulse_response.*hanning(max(size(xmt_impulse_response)));
rcv_impulse_response = xmt_impulse_response;

excitation = (sin(2*pi*f0*(0:1/fs:2/f0)))';
excitation = excitation.*hanning(max(size(excitation)));
                                    
%% Generate scatter field

sca_x = ([ 0  0  0  0  0  0  0  0  0]./1000)';
sca_y = ([zeros(size(sca_x,1),1)]./1000);
sca_z = ([15 25 35 45 55 65 75 85 95]./1000)';

media.phantom_positions = [sca_x sca_y sca_z];
media.phantom_amplitudes = ones(size(media.phantom_positions,1),1);


%%

%% Generate Data
create_data = 1;
if(create_data == 1)    
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

end
      return
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define second stage beamforming parameters And view parameters
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

useCaseParams.bfrcvparams(1).rcvfnum = useCaseParams.bfxmitparams(1).xmitfnum;
useCaseParams.bfrcvparams(1).rcvapodishape = 1;
useCaseParams.bfrcvparams(1).rcvapodigausswidth = 0.2;
useCaseParams.bfrcvparams(1).rcvfocus = useCaseParams.bfxmitparams(1).xmitfocus;

useCaseParams.scanparams(1).startdepthq = 0.0;
useCaseParams.scanparams(1).stopdepthq = 0.12;

useCaseParams.scanparams(1).startlinenumq = 0;
useCaseParams.scanparams(1).stoplinenumq = 268;

useCaseParams.scanparams(1).scanareadscr.startlineorigin.y = 0;
useCaseParams.scanparams(1).scanareadscr.startlineorigin.x = 0;
useCaseParams.scanparams(1).scanareadscr.startlineangle = pi/2-30/180*pi;

useCaseParams.scanparams(1).scanareadscr.stoplineorigin.y = 0;
useCaseParams.scanparams(1).scanareadscr.stoplineorigin.x = 0;
useCaseParams.scanparams(1).scanareadscr.stoplineangle =  pi/2+30/180*pi;

useCaseParams.scanparams(1).steeringangle = 0;
useCaseParams.scanparams(1).phasedangle = 30/180*pi;

% View parameters - included here because they must be saved with usecase
useCaseParams.scanparams(1).windowtissueq.x_tismin = -0.04;
useCaseParams.scanparams(1).windowtissueq.y_tismin = 0;                
useCaseParams.scanparams(1).windowtissueq.x_tismax = 0.04;
useCaseParams.scanparams(1).windowtissueq.y_tismax = 0.12;
useCaseParams.scanparams(1).scantype = 2;
useCaseParams.bfrcvparams(1).smpfreq = 30e6;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define second stage beamforming parameters And view parameters
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%% Beamform 2nd stage
tic
BeamformationVer4('beamformation_method','sasbsim_new', ...
                  'usecase',useCaseParams, ...
                  'loadpath', savepath, ...
                  'symmetric','symmetric',...
                  'savepath',[savepath],...
                  'generate_scaling',false)
toc           

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot
%% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_nr = 8;
type = 'B-mode';
switch(type)
    case{'psf'}
        compression = 'Linear';
        dynamic_range = 60;
    case{'B-mode'}
        compression = 'BK_muLaw_RTSC_VP_DRC_SEL_0';
        dynamic_range = 100;
end
plot_obj_sasb_2 = Plot_Beamformed_Data_Ver2(...
        'figure nr', 3, ...
        'data type','simulation',...
        'type of beamformation', 'sasb',...
        'stage','second',...
        'show tgc','dont show tgc',...
        'loadpath',[savepath ],...
        'gain',[],...
        'windowtissueq',useCaseParams.scanparams(1).windowtissueq,...
        'compression',compression,...
        'dynamic range',dynamic_range); 
    
    
set(gcf,'position',[1922         506         560         420])
drawnow
axis image
caxis([-60 0])




