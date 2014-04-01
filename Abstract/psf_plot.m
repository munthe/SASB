% Script to plot psf measurements for use in abstract for ieee ultrasound
% symposium

%% Initilize
base_path = './';
loadpath = './';
base_path_usecase = './';
addpath(genpath('../lib'))
savepath = loadpath;
savepath_figures = savepath;
figure_format = '-depsc2';
setupfigdefault

%% Load point spread function data, output from cfu_get_psf_metrics
load psf

%% Plot psf as function of depth
n = 1:7;

fig_nr = 1;
set(0, 'DefaultAxesColorOrder', linspecer(2));
figure(fig_nr);
depth = cell2mat({psf_SMF(n).y_coord});
plot( ...
    depth,cell2mat({psf_SMF(n).fwhm_x}),'^-', ...
    depth,cell2mat({psf_SASB(n).fwhm_x}),'v-', ...
    depth,cell2mat({psf_SMF(n).radius20dB}),'^--', ...
    depth,cell2mat({psf_SASB(n).radius20dB}),'v--',...
    'MarkerSize',6, ...
    'LineWidth',1.5 ...
    );
xlabel('Depth of point scatter [mm]');
ylabel('[mm]');
legend('FWHM SMF','FWHM SASB','R\_20dB SMF','R\_20dB SASB','Location','NorthWest');
% title(psf(setting,1).setupDesc);
prettyfig
print(figure(fig_nr),[savepath_figures 'psf'],figure_format)
hold off

%% Calculate percentages

mean(cell2mat({psf_SMF(n).fwhm_x})./cell2mat({psf_SASB(n).fwhm_x}))
mean(cell2mat({psf_SMF(n).fwhm_x}))/mean(cell2mat({psf_SASB(n).fwhm_x}))

mean(cell2mat({psf_SMF(n).radius20dB})./cell2mat({psf_SASB(n).radius20dB}))
mean(cell2mat({psf_SMF(n).radius20dB}))/mean(cell2mat({psf_SASB(n).radius20dB}))