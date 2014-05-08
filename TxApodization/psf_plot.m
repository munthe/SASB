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

%% Plot psf as function of depth, fwhm
n = 1:9;
c = linspecer(5);
fig_nr = 1;
set(0, 'DefaultAxesColorOrder', [c(2,:);c(1,:);c(5,:);c(3,:)]);
figure(fig_nr);
depth = cell2mat({psf(1,:).y_coord});
plot1 = plot( ...
    depth,cell2mat({psf_SMF.fwhm_x}),'^-', ...
    depth,cell2mat({psf(1,:).fwhm_x}),'v-', ...
    depth,cell2mat({psf(14,:).fwhm_x}),'o-', ...
    depth,cell2mat({psf_SMFTx(1,:).fwhm_x}),'*-', ...
    'MarkerSize',6, ...
    'LineWidth',1 ...
    );
set(plot1(1),...
     'MarkerFaceColor',c(2,:));
set(plot1(2),...
     'MarkerFaceColor',c(1,:));
set(plot1(3),...
     'MarkerFaceColor',c(5,:));
set(plot1(4),...
     'MarkerFaceColor',c(3,:));

xlabel('Depth of point scatter [mm]');
ylabel('[mm]');
axis([10 100 1.2 2.5]);

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
legend('SMF', 'SASB', '0,0.6,1','SMF+TxApo','Location','NorthWest');
legend boxoff
% title(psf(setting,1).setupDesc);
prettyfig(14)

print(figure(fig_nr),[savepath_figures 'psf_FWHM'],figure_format)
hold off

%% Plot psf as function of depth, radius20dB
n = 1:9;
c = linspecer(5);
fig_nr = 2;
set(0, 'DefaultAxesColorOrder', [c(2,:);c(1,:);c(5,:);c(3,:)]);
figure(fig_nr);
depth = cell2mat({psf(1,:).y_coord});
plot1 = plot( ...
    depth,cell2mat({psf_SMF.radius20dB}),'^--', ...
    depth,cell2mat({psf(1,:).radius20dB}),'v--',...
    depth,cell2mat({psf(14,:).radius20dB}),'o--',...
    depth,cell2mat({psf_SMFTx(1,:).radius20dB}),'*--', ...
    'MarkerSize',6, ...
    'LineWidth',1 ...
    );

set(plot1(1),...
      'MarkerFaceColor',c(2,:));
set(plot1(2),...
     'MarkerFaceColor',c(1,:));
set(plot1(3),...
     'MarkerFaceColor',c(5,:));
set(plot1(4),...
     'MarkerFaceColor',c(3,:));

xlabel('Depth of point scatter [mm]');
ylabel('[mm]');
axis([10 100 1.4 3]);

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
legend('SMF', 'SASB', '0,0.6,1','SMF+TxApo','Location','NorthWest');
legend boxoff
% title(psf(setting,1).setupDesc);
prettyfig(14)

print(figure(fig_nr),[savepath_figures 'psf_R20dB'],figure_format)
hold off

%% Plot cystic resolution for scatter 7
n = 7;
c = linspecer(5);
fig_nr = 3;
set(0, 'DefaultAxesColorOrder', [c(2,:);c(1,:);c(5,:);c(3,:)]);

figure(fig_nr)
plot(...
    psf_SMF(1,n).radius,psf_SMF(1,n).ct,...
    psf(1,n).radius,psf(1,n).ct,...
    psf(14,n).radius,psf(14,n).ct,...
    psf_SMFTx(1,n).radius,psf_SMFTx(1,n).ct,...
    'LineWidth',1.5 ....
    );
xlabel('Radius [mm]')
ylabel('Attenuation [dB]')

legend('SMF', 'SASB, boxcar', '0,0.6,1','SMF+TxApo','Location','NorthEast');
legend boxoff
axis([0 9.5 -35 0]);

prettyfig(10)
print(figure(fig_nr),[savepath_figures 'cysticresolution'],figure_format)

%% Calculate percentages

fprintf('FWHM reduced with %f',...
    1-mean(cell2mat({psf_SMF(1,:).fwhm_x})./cell2mat({psf(14,:).fwhm_x})))
    %mean(cell2mat({psf(14,:).fwhm_x}))/mean(cell2mat({psf(1,:).fwhm_x})))

fprintf('R_{20dB} increased to %f',...
    mean(cell2mat({psf_SMF(1,:).radius20dB})./cell2mat({psf(14,:).radius20dB})))
    %1-mean(cell2mat({psf(14,:).radius20dB}))/mean(cell2mat({psf(1,:).radius20dB})))

