function handles = plot_apo (obj, fontsize)
% PLOT_APO
% 
% This function overlays the apodization values of each element on the transducer.
% 
% It accepts one input argument: the font size.
% It outputs an optional struct: the colorbar and figure axis handles.
%
%  ems_obj = mfr_emission('c',1504, 'no_act_elm',512);
%  fontsize = 16;
%  figure;
%  handles = ems_obj.plot_apo(fontsize);
%  set(handles.fig_ax, 'YTick', [])
%  ylabel('')
%
%  'fontsize' defaults to 20.
%
% By MOFI, 2011, 2012, 2013.

if nargin < 2
    fontsize = 20;
end


% Get element limits
x_min = min(obj.pos(:,1))-obj.pitch_x/2;
x_max = max(obj.pos(:,1))+obj.pitch_x/2;
y_min = min(obj.pos(:,2))-obj.pitch_y/2;
y_max = max(obj.pos(:,2))+obj.pitch_y/2;


% make background that will also export
bck_x = [x_min-obj.pitch_x x_max+obj.pitch_x];
bck_y = [y_min-obj.pitch_y y_max+obj.pitch_y];
[x y z] = meshgrid(bck_x*1e3,bck_y*1e3,0);
c = 0;
%surf(x,y,z,c);
%hold on

% Plot a surf of each element
for idx = 1:obj.no_elm
    pos = obj.pos(idx,:);
    elm_x = [pos(1)-obj.pitch_x/2 pos(1)+obj.pitch_x/2];
    elm_y = [pos(2)-obj.pitch_y/2 pos(2)+obj.pitch_y/2];
    [x y z] = meshgrid(elm_x*1000,elm_y*1000,0);
    c = obj.apo(idx);
    %c = c *(0.995)+0.005;
    surf(x,y,z,c);
    hold on;
end
%plot3(obj.origin_coord(1)*1000, obj.origin_coord(2)*1000, 0.1*1000, 'k*', 'MarkerSize',12)
hold off;


% Set axes and labels
view (0,90)
axis equal;
xlim([(x_min-obj.pitch_x)*1e3 (x_max+obj.pitch_x)*1e3]);
ylim([(y_min-obj.pitch_y)*1e3 (y_max+obj.pitch_y)*1e3]);
ca = gca; %axis handle
set(ca,'Xtick',[x_min*1e3 0 x_max*1e3])
set(ca,'Ytick',[y_min*1e3 0 y_max*1e3])
xlabel('x-axis [mm]', 'FontSize',fontsize)
ylh = ylabel('y-axis [mm]', 'FontSize',fontsize);
% $$$ ylpos = get(ylh,'Position');
% $$$ ylpos = ylpos + [1.5 0 0];
% $$$ set(ylh,'Position', ylpos)

set(ca, 'FontSize',fontsize)


% Colorbar + Colormap
caxis([0 1])
cmap = flipud(thermal(200));
% $$$ cmap(2,:) = [1 1 1];
% $$$ cmap(1,:) = [0.8 0.8 0.8];
colormap (cmap)
cbh = colorbar;
set(cbh, 'Ytick',[0 0.98]) %lower the '1' a bit.
set(cbh, 'YtickLabel', ['0'; '1'])

set(cbh, 'FontSize',fontsize)
ylabel(cbh, 'Apodization lvl.', 'FontSize',fontsize)
ylpos = get(get(cbh,'Ylabel'),'Position');
ylpos = ylpos - [2 0 0];
set(get(cbh,'Ylabel'),'Position', ylpos)

% set background color
set(gcf, 'InvertHardCopy', 'off'); % keeps background color when exporting
set(gca,'Color', [0.8 0.8 0.8])
set(gcf,'Color', [1 1 1])


grid off;
box on;

% set output
if nargout > 0
    handles.cb_ax  = cbh;
    handles.fig_ax = ca;
end

