function handles = plot_vs (obj)
% PLOT_VS
% 
% This function plots the vs and overlays the apodization values of each element on the transducer.
% 
% It accepts one input argument: the font size.
% It outputs an optional struct containing the colorbar and figure axis handles.
%
%  ems_obj = mfr_emission('c',1504, 'no_act_elm',512);
%  fontsize = 16;
%  figure;
%  handles = ems_obj.plot_vs(fontsize);
%  set(handles.fig_ax, 'YTick', [])
%  ylabel('')
%
%  'fontsize' defaults to 20.
%
% By MOFI, 2011, 2012, 2013.

if nargin < 2
    fontsize = 20;
end


pos = obj.pos(1,:);
elm_x = [pos(1)-obj.pitch_x/2 pos(1)+obj.pitch_x/2];
elm_y = [pos(2)-obj.pitch_y/2 pos(2)+obj.pitch_y/2];
%elm_z = [0 0];
% $$$ c = obj.delays(1)*obj.c;
% $$$ [x y z] = meshgrid(elm_x*1000,elm_y*1000,-c*1000);
c = obj.apo(1);
[x y z] = meshgrid(elm_x*1000,elm_y*1000,0);
surf(x,y,z,c);
% $$$ if extra_opt
% $$$     surf(x,y,z,c, plot_opt_var, plot_opt_val);
% $$$ else
% $$$     surf(x,y,z,c);
% $$$ end

hold on;

for idx = 2:obj.no_elm
    pos = obj.pos(idx,:);
    elm_x = [pos(1)-obj.pitch_x/2 pos(1)+obj.pitch_x/2];
    elm_y = [pos(2)-obj.pitch_y/2 pos(2)+obj.pitch_y/2];
    %elm_z = [0 0];
    %c = obj.delays(idx)*obj.c;
    %[x y z] = meshgrid(elm_x*1000,elm_y*1000,-c*1000);
    c = obj.apo(idx);
    [x y z] = meshgrid(elm_x*1000,elm_y*1000,0);
    surf(x,y,z,c);
% $$$     if extra_opt
% $$$         surf(x,y,z,c, plot_opt_var, plot_opt_val);
% $$$     else
% $$$     end
end
 
plot3(obj.vs(1)*1000, obj.vs(2)*1000, obj.vs(3)*1000, '*')
plot3(obj.origin_coord(1)*1000, obj.origin_coord(2)*1000, obj.origin_coord(3)*1000, 'g*')
% $$$ if extra_opt
% $$$     plot3(obj.vs(1)*1000, obj.vs(2)*1000, obj.vs(3)*1000, '*', plot_opt_var, plot_opt_val)
% $$$ else
% $$$     plot3(obj.vs(1)*1000, obj.vs(2)*1000, obj.vs(3)*1000, '*')
% $$$     plot3(obj.origin_coord(1)*1000, obj.origin_coord(2)*1000, obj.origin_coord(3)*1000, 'g*')
% $$$ end

plot3([obj.origin_coord(1) obj.vs(1)]*1e3, ...
      [obj.origin_coord(1) obj.vs(2)]*1e3, ...
      [obj.origin_coord(1) obj.vs(3)]*1e3, 'b-')
    
    
axis equal;
%view (0,90)
xlabel('x-axis [mm]')
ylabel('y-axis [mm]')
zlabel('z-axis [mm]')


hold off;

% Get element limits
x_min = min([obj.pos(:,1); obj.vs(1)])-obj.pitch_x/2;
x_max = max([obj.pos(:,1); obj.vs(1)])+obj.pitch_x/2;
y_min = min([obj.pos(:,2); obj.vs(2)])-obj.pitch_y/2;
y_max = max([obj.pos(:,2); obj.vs(2)])+obj.pitch_y/2;
% $$$ z_min = min([0 obj.vs(3)]*1e3);
% $$$ z_max = max([0 obj.vs(3)]*1e3);
z_min = min([0 obj.focus_r]);
z_max = max([0 obj.focus_r]);


% Set axes and labels
%view (0,90)
axis equal;
xlim([(x_min-obj.pitch_x)*1e3 (x_max+obj.pitch_x)*1e3]);
ylim([(y_min-obj.pitch_y)*1e3 (y_max+obj.pitch_y)*1e3]);
zlim([z_min z_max]*1e3);

ca = gca; %axis handle
set(ca,'Xtick',[x_min*1e3 0 x_max*1e3])
set(ca,'Ytick',[y_min*1e3 0 y_max*1e3])

xlabel('x-axis [mm]', 'FontSize',fontsize)
ylabel('y-axis [mm]', 'FontSize',fontsize)
set(ca, 'FontSize',fontsize)



% Colorbar + Colormap
caxis([0 1])
cmap = flipud(thermal(200));
cmap(1,:) = [0.95 0.95 0.95];
colormap (cmap)
% $$$ cbh = colorbar;
% $$$ set(cbh, 'Ytick',[0 0.96])
% $$$ set(cbh, 'YtickLabel', ['0'; '1'])
% $$$ 
% $$$ set(cbh, 'FontSize',fontsize)
% $$$ ylabel(cbh, 'Apodization lvl.', 'FontSize',fontsize)
% $$$ ylpos = get(get(cbh,'Ylabel'),'Position');
% $$$ ylpos = ylpos - [2 0 0];
% $$$ set(get(cbh,'Ylabel'),'Position', ylpos)

grid off;
box on;

% set output
if nargout > 0
    handles.cb_ax  = cbh;
    handles.fig_ax = ca;
end

