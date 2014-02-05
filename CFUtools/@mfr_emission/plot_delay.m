function handles = plot_delay (obj, fontsize)
% PLOT_APO
% 
% This function overlays the delay value of each element on the transducer.
% 
% It accepts one input argument: the font size.
% It outputs an optional struct: the colorbar and figure axis handles.
%
%  ems_obj = mfr_emission('c',1504, 'no_act_elm',512);
%  fontsize = 16;
%  figure;
%  handles = ems_obj.plot_delay(fontsize);
%  set(handles.fig_ax, 'YTick', [])
%  ylabel('')
%
%  'fontsize' defaults to 20.
%
% By MOFI, 2011, 2012, 2013.

if nargin < 2
    fontsize = 20;
end

% Plot a surf of each element
for idx = 1:obj.no_elm
    pos = obj.pos(idx,:);
    elm_x = [pos(1)-obj.pitch_x/2 pos(1)+obj.pitch_x/2];
    elm_y = [pos(2)-obj.pitch_y/2 pos(2)+obj.pitch_y/2];
    [x y z] = meshgrid(elm_x*1000,elm_y*1000,0);
    c = obj.delays(idx)*1e6;
    surf(x,y,z,c);
    hold on;
end
%plot3(obj.origin_coord(1)*1000, obj.origin_coord(2)*1000, 0.1*1000, 'k*', 'MarkerSize',12)
hold off;


% Get element limits
x_min = min(obj.pos(:,1))-obj.pitch_x/2;
x_max = max(obj.pos(:,1))+obj.pitch_x/2;
y_min = min(obj.pos(:,2))-obj.pitch_y/2;
y_max = max(obj.pos(:,2))+obj.pitch_y/2;


% Set axes and labels
view (0,90)
axis equal;
xlim([(x_min-obj.pitch_x)*1e3 (x_max+obj.pitch_x)*1e3]);
ylim([(y_min-obj.pitch_y)*1e3 (y_max+obj.pitch_y)*1e3]);
ca = gca; %axis handle
set(ca,'Xtick',[x_min*1e3 0 x_max*1e3])
set(ca,'Ytick',[y_min*1e3 0 y_max*1e3])
xlabel('x-axis [mm]', 'FontSize',fontsize)
ylabel('y-axis [mm]', 'FontSize',fontsize)
set(ca, 'FontSize',fontsize)


% Colorbar + Colormap
cbh = colorbar;
max_delay = max(obj.delays*1e6);
cmap = flipud(thermal(200));
%cmap(1,:) = [0.95 0.95 0.95];
colormap (cmap)
set(cbh, 'FontSize',fontsize)
caxis([0 max_delay])
max_delay = floor(max_delay*100)/100; %max two decimals
set(cbh, 'Ytick',[0 max_delay])
ylabel(cbh, 'Delay [µs]', 'FontSize',fontsize)
ylpos = get(get(cbh,'Ylabel'),'Position');
ylpos = ylpos - [5 0 0];
set(get(cbh,'Ylabel'),'Position', ylpos)

grid off;
box on;
set(gca, 'Color',[0.8 0.8 0.8])
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'Color', [1 1 1]);


% set output
if nargout > 0
    handles.cb_ax  = cbh;
    handles.fig_ax = ca;
end
    
