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
    [-1,-2],[1,1],'-',...
    depth,cell2mat({psf_SMF(n).fwhm_x}),'^-', ...
    depth,cell2mat({psf_SASB(n).fwhm_x}),'v-', ...
    depth,cell2mat({psf_SMF(n).radius20dB}),'^--', ...
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
axis([10 80 1.2 2.7]);

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
legend('SMF', 'SASB','Location','NorthWest');
legend boxoff
% title(psf(setting,1).setupDesc);
prettyfig
print(figure(fig_nr),[savepath_figures 'psf'],figure_format)
hold off

%% Calculate percentages

fprintf('FWHM reduced with %f or %f\n',...
    1-mean(cell2mat({psf_SMF(n).fwhm_x})./cell2mat({psf_SASB(n).fwhm_x})),...
    1-mean(cell2mat({psf_SMF(n).fwhm_x}))/mean(cell2mat({psf_SASB(n).fwhm_x})))

fprintf('R_{20dB} increased to %f or %f\n',...
    mean(cell2mat({psf_SMF(n).radius20dB})./cell2mat({psf_SASB(n).radius20dB})),...
    mean(cell2mat({psf_SMF(n).radius20dB}))/mean(cell2mat({psf_SASB(n).radius20dB})))