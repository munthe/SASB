% This file gives an example of how to use the cfu_video*
% functions.
% By MOFI
% V. 1.0, tor 02 maj 2013 00:07:26 CEST
% V. 1.1, fre 03 maj 2013 11:21:25 CEST, changed frame rate to 20.53, which results 
%                                        in a frame rate of 23.97.






%% Setup cfu_video
% Only 'video_name' is mandatory to set. See the help file for
% default values.
video_name     = 'double_circle2';
video_settings = cfu_video_init(video_name, ...
                                'frame_rate', 20.53, ...
                                'folder','video_test_mpeg1',...
                                'visible', 'off', ...
                                'author', 'Morten Fischer Rasmussen', ...
                                'resolution', [1280 960], ...
                                'encoding', 'MPEG1');


%% Setup plot
circle1_radius = 3;
circle2_radius = 1;

t1 = (0:0.05:2*pi)';
t2 = (0:0.25:10*pi)';

circle1 = [cos(t1) sin(t1)]*circle1_radius;
circle2 = [cos(t2) sin(t2)]*circle2_radius;

wobble = circle1 + circle2;

fontsize = 25;
% Increase the axes font size
set(video_settings.hax,'FontSize', fontsize)


%% Plot + grab frame
for idx = 1:length(t1)
    fprintf('frame: %i/%i\n', idx, length(t1));
    % clear figure
    % setup plot coordinates
    pt_x   = wobble(idx,1);
    pt_y   = wobble(idx,2);
    line_x = wobble(1:idx,1);
    line_y = wobble(1:idx,2);
    
    %plot
    hp1 = plot(video_settings.hax, line_x,line_y, 'k-');
    hold on
    hp2 = plot(video_settings.hax, pt_x,pt_y, 'r.', 'MarkerSize',20 );
    axis image
    axis tight
    xlim([-4.1 4.1])
    ylim([-4.1 4.1])
    xlabel('x [mm]', 'FontSize', fontsize)
    ylabel('y [mm]', 'FontSize', fontsize)
    grid on;
    box on;
    drawnow;
    
    cfu_video_grab_frame(video_settings, idx);
    % Delete the plot without deleting the axes (don't use clf!)
    delete(hp1)
    delete(hp2)
end





%% encode video
% Delete the temporary images from the output folder
delete_tmp_images = 0;
cfu_video_encode(video_settings, delete_tmp_images);

