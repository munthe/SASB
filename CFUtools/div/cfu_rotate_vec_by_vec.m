function v_rot = cfu_rotate_vec_by_vec(v, rot_vec)
% v_rot = cfu_rotate_vec_by_vec(v, rot_vec)
%
% Rotates the [3x1] vector 'v' by the unit vector 'rot_vec'. 
% 'v' is rotated in the same direction and by the same amount as 'rot_vec' is rotated 
% away from the z-axis ([0;0;1]).
% 
% By MOFI
% 2013-12-19, Version 1.0, Init version.
%
%

% make sure the rotation vector is a unit vector
rot_vec = rot_vec./norm(rot_vec); 

z_axis = [0;0;1];

% angle between the rotation vector ('rot_vec') and the z-axis
phi = -acos(dot(rot_vec, z_axis))

% axis of rotation
rot_axis = cross(rot_vec, z_axis)
if norm(rot_axis) ==0 % vector is parallel to z_axis
    rot_axis = [1;0;0];
end
% normalise the vector
L0 = rot_axis/norm(rot_axis);

% Rodrigues' rotation formula, from http://en.wikipedia.org/wiki/Axis-angle:
v_rot = v*cos(phi) + ...
        (v'*L0)*(1-cos(phi))*L0 + ...
        cross(L0,v)*sin(phi);

