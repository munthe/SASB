%  Create a computer model of a point phantom. The phantom contains
%  20 point targets separated by 5 mm and starting 10 mm from the
%  transducer surface.
%
%  Calling: [positions, amp] = pts_pha ;
%
%  Output:      positions  - Positions of the scatterers.
%               amp        - amplitude of the scatterers.
%
%  Version 1.0, March 25, 1997 by Joergen Arendt Jensen

function [positions, amp] = point_phantom 

dz=5/1000;        %  Distance between points [m]
z_start=10/1000;  %  Start of point [m]
Npoints=20;       %  Number of point targets

%  Create the point scatterers positions

positions = [zeros(1,Npoints); zeros(1,Npoints); (1:Npoints)*dz+z_start]';
amp=ones(Npoints,1);

end
