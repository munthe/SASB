% Script to plot psf measurements for use in abstract for ieee ultrasound
% symposium

%% Initilize
base_path = './';
loadpath = './';
base_path_usecase = './';

savepath = loadpath;
savepath_figures = savepath;
figure_format = '-depsc2';
setupfigdefault

%% Load point spread function data, output from cfu_get_psf_metrics
load psf

%% Plot psf as function of depth

fig_nr = 1;
figure(fig_nr);
n = 1:7;

depth = cell2mat({psf_SMF(n).y_coord})*1000;
plot( ...
    depth,cell2mat({psf_SMF(n).fwhm_x}),'o--', ...
    depth,cell2mat({psf_SASB(n).fwhm_x}),'d--', ...
    depth,cell2mat({psf_SMF(n).radius20dB}),'o:', ...
    depth,cell2mat({psf_SASB(n).radius20dB}),'d:' ...
    );
xlabel('Depth of point scatter [mm]');
ylabel('[mm]');
legend('FWHM SMF','FWHM SASB','Radius 20dB SMF','Radius 20dB SASB','Location','NorthWest');
% title(psf(setting,1).setupDesc);
prettyfig
print(figure(fig_nr),[savepath_figures 'psf'],figure_format)