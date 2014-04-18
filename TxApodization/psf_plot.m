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
load Transducer8804_Focus0.02_F#2.mat

%% Plot psf as function of depth
n = 1:9;
c = linspecer(5);
fig_nr = 1;
set(0, 'DefaultAxesColorOrder', [c(1,:);c(2,:);c(5,:)]);
figure(fig_nr);
depth = cell2mat({psf(1,:).y_coord});
plot1 = plot( ...
    depth,cell2mat({psf(1,:).fwhm_x}),'^-', ...
    depth,cell2mat({psf(2,:).fwhm_x}),'v-', ...
    depth,cell2mat({psf(14,:).fwhm_x}),'o-', ...
    depth,cell2mat({psf(1,:).radius20dB}),'^--', ...
    depth,cell2mat({psf(2,:).radius20dB}),'v--',...
    depth,cell2mat({psf(14,:).radius20dB}),'o--',...
    'MarkerSize',6, ...
    'LineWidth',1 ...
    );
set(plot1(1),...
     'MarkerFaceColor',c(1,:));
set(plot1(2),...
     'MarkerFaceColor',c(2,:));
set(plot1(3),...
     'MarkerFaceColor',c(5,:));
set(plot1(4),...
      'MarkerFaceColor',c(1,:));
set(plot1(5),...
     'MarkerFaceColor',c(2,:));
set(plot1(6),...
     'MarkerFaceColor',c(5,:));

xlabel('Depth of point scatter [mm]');
ylabel('[mm]');
axis([10 100 1.4 2.6]);

% text(depth(7),(cell2mat({psf_SMF(7).fwhm_x})+cell2mat({psf_SMF(7).radius20dB}))/2,'SMF','HorizontalAlignment','left');
% %text(depth(1),cell2mat({psf_SMF(7).radius20dB})+0.1,'SMF','HorizontalAlignment','center');
% text(depth(7),(cell2mat({psf_SASB(7).fwhm_x})+cell2mat({psf_SASB(7).radius20dB}))/2,'SASB','HorizontalAlignment','left');
% %text(depth(1),cell2mat({psf_SASB(7).radius20dB})+0.1,'SASB','HorizontalAlignment','center');

% Create arrow
% annotation('arrow',[0.84 0.81],...
%     [0.783444444444444 0.855555555555556],'Color',[0.3 0.3 0.3]);
% annotation('arrow',[0.84 0.81],...
%     [0.359 0.222222222222222],'Color',[0.3 0.3 0.3]);
% annotation('arrow',[0.84 0.81],...
%     [0.425666666666667 0.526666666666667],'Color',[0.3 0.3 0.3]);
% annotation('arrow',[0.84 0.81],...
%     [0.725666666666667 0.628888888888889],'Color',[0.3 0.3 0.3]);

% legend('FWHM', 'R_{20dB}','Location','NorthWest');
legend('Boxcar', 'Inf', '0,0.6,1','Location','NorthWest');
legend boxoff
% title(psf(setting,1).setupDesc);
prettyfig(8)
print(figure(fig_nr),[savepath_figures 'psf'],figure_format)
hold off

%% Calculate percentages

fprintf('FWHM reduced with %f or %f\n',...
    1-mean(cell2mat({psf_SMF(n).fwhm_x})./cell2mat({psf_SASB(n).fwhm_x})),...
    1-mean(cell2mat({psf_SMF(n).fwhm_x}))/mean(cell2mat({psf_SASB(n).fwhm_x})))

fprintf('R_{20dB} increased to %f or %f\n',...
    mean(cell2mat({psf_SMF(n).radius20dB})./cell2mat({psf_SASB(n).radius20dB})),...
    mean(cell2mat({psf_SMF(n).radius20dB}))/mean(cell2mat({psf_SASB(n).radius20dB})))