function [rect] = Create_Linear_Focused_Element(obj,center,pitch,kerf,height,elevation_focus,aclayerthickness,nr_sub_x,nr_sub_y)

% calculate width of physical element
width = pitch-kerf;


% determine points in azimut 
x = linspace(center(1)-width/2,center(1)+width/2,nr_sub_x+1)';
z = repmat(center(3),size(x,1),1);
y = zeros(size(x,1),1);
pos_top = [x y z];

% determine points for elevation arc
StartAngle_inner = (acos(height/2/elevation_focus));
StopAngle_inner = pi-StartAngle_inner;
[y_inner_arc,z_inner_arc] = Calculate_arc_points(elevation_focus,elevation_focus,StartAngle_inner,StopAngle_inner,nr_sub_y+1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Added 18/3-2011 by Martin C. Hemmsen
% pos_top(:,3) = pos_top(:,3)+z_inner_arc(1); 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% StartAngle_inner = -(pi-pi/2-acos(height/2/elevation_focus));
% StopAngle_inner = -StartAngle_inner;
% [y_inner_arc,z_inner_arc] = Calculate_arc_points(elevation_focus,elevation_focus,StartAngle_inner,StopAngle_inner,nr_sub_y+1);
% 
% figure
% plot3(pos_top(:,1),pos_top(:,2),pos_top(:,3),'-or')
% hold on

% calculate corners etc.
rect = zeros(nr_sub_y*nr_sub_x,19);
for k = 1:nr_sub_x 
    pos_right_arc = [repmat(pos_top(k,1),nr_sub_y+1,1) y_inner_arc' ...
        pos_top(k,3)-z_inner_arc'];
    pos_left_arc = [repmat(pos_top(k+1,1),nr_sub_y+1,1) y_inner_arc' ...
        pos_top(k+1,3)-z_inner_arc'];
%     
%     % for plotting
%     plot3(pos_left_arc(:,1),pos_left_arc(:,2),pos_left_arc(:,3))
%     plot3(pos_right_arc(:,1),pos_right_arc(:,2),pos_right_arc(:,3),'r')
%     % end for plotting
    
    % build corners
    corners = zeros(nr_sub_y,3*4);
    mat_element_center = zeros(nr_sub_y,3);
    for n = 1:nr_sub_y
        corners(n,1:3) = pos_left_arc(n,:);
        corners(n,4:6) = pos_left_arc(n+1,:);
        corners(n,7:9) = pos_right_arc(n+1,:);
        corners(n,10:12) = pos_right_arc(n,:);
        
        
        mat_element_center(n,:) = (pos_left_arc(n,:)-pos_right_arc(n+1,:))/2+pos_right_arc(n+1,:);
        
%         % for plotting
%         vert  = [corners(n,1:3);
%                  corners(n,4:6);
%                  corners(n,7:9);
%                  corners(n,10:12)];
%         fac = [1 2 3; 1 3 4];
%         patch('Faces',fac,'Vertices',vert)
%         
%         plot3(mat_element_center(n,1),mat_element_center(n,2),mat_element_center(n,3),'ob')
%         plot3(corners(n,1),corners(n,2),corners(n,3),'xr','linewidth',2)
%         plot3(corners(n,4),corners(n,5),corners(n,6),'xb','linewidth',2)
%         plot3(corners(n,7),corners(n,8),corners(n,9),'xk','linewidth',2)
%         plot3(corners(n,10),corners(n,11),corners(n,12),'xg','linewidth',2)
%         % end for plotting
    end
    
    % build rects
    rect(1+(k-1)*nr_sub_y:k*nr_sub_y,1) = 0;  % The number for the physical aperture starting from one
    rect(1+(k-1)*nr_sub_y:k*nr_sub_y,2:13) = corners;
    rect(1+(k-1)*nr_sub_y:k*nr_sub_y,14) = 1; % Apodization value for this element.
    rect(1+(k-1)*nr_sub_y:k*nr_sub_y,15) = width; % Width of the element (x direction)    
    rect(1+(k-1)*nr_sub_y:k*nr_sub_y,16) = height; % Height of the element (y direction)
    rect(1+(k-1)*nr_sub_y:k*nr_sub_y,17:19) = mat_element_center;

    
end

% % for plotting
% xlabel('x')
% ylabel('y')
% zlabel('z')
% 
% 
% center=[0 0 0];
% focus=[0 0 70]/1000;
% Th = xdc_rectangles (rect, center, focus);
% % %xdc_show(Th)
% figure
% show_xdc(Th)
% xdc_free(Th)
% % end for plotting



