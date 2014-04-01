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
c = linspecer(2);
fig_nr = 1;
set(0, 'DefaultAxesColorOrder', c);
figure(fig_nr);
depth = cell2mat({psf_SMF(n).y_coord});
plot1 = plot( ...
    [-1,-2],[0,0],'-',...
    [-1,-2],[1,1],'--',...
    depth,cell2mat({psf_SMF(n).fwhm_x}),'^-', ...
    depth,cell2mat({psf_SMF(n).radius20dB}),'^--', ...
    depth,cell2mat({psf_SASB(n).fwhm_x}),'v-', ...
    depth,cell2mat({psf_SASB(n).radius20dB}),'v--',...
    'MarkerSize',6, ...
    'LineWidth',1.5 ...
    );
set(plot1(3),...
     'MarkerFaceColor',c(1,:));
 set(plot1(4),...
     'MarkerFaceColor',c(2,:));
 set(plot1(5),...
      'MarkerFaceColor',c(1,:));
  set(plot1(6),...
     'MarkerFaceColor',c(2,:));
xlabel('Depth of point scatter [mm]');
ylabel('[mm]');
axis([10 80 0.5 2.5]);
text(depth(1),cell2mat({psf_SMF(1).fwhm_x})-0.1,'SMF','HorizontalAlignment','center');
%text(depth(1),cell2mat({psf_SMF(1).radius20dB})+0.1,'SMF','HorizontalAlignment','center');
text(depth(1),cell2mat({psf_SASB(1).fwhm_x})-0.1,'SASB','HorizontalAlignment','center');
%text(depth(1),cell2mat({psf_SASB(1).radius20dB})+0.1,'SASB','HorizontalAlignment','center');

legend('FWHM', 'R\_20dB','Location','NorthWest');
% title(psf(setting,1).setupDesc);
prettyfig
print(figure(fig_nr),[savepath_figures 'psf'],figure_format)
hold off

%% Calculate percentages

mean(cell2mat({psf_SMF(n).fwhm_x})./cell2mat({psf_SASB(n).fwhm_x}))
mean(cell2mat({psf_SMF(n).fwhm_x}))/mean(cell2mat({psf_SASB(n).fwhm_x}))

mean(cell2mat({psf_SMF(n).radius20dB})./cell2mat({psf_SASB(n).radius20dB}))
mean(cell2mat({psf_SMF(n).radius20dB}))/mean(cell2mat({psf_SASB(n).radius20dB}))