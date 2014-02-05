function [euler basis cart] = cfu_beam2euler_angles (zx, zy)
%
% euler = cfu_beam2euler_angles (zy, zx)
% euler = [theta psi phi]
%
% 2011-11-20
% MFR


% get cartesian coordinate
cart = cfu_beam2cart_coord(1,zy,zx);

% when only the two first rotations are used, euler and spherical angles are equal
% see: http://en.wikipedia.org/wiki/File:EulerProjections.svg or
%      http://en.wikipedia.org/wiki/File:EulerG.png

% cartesian to spherical 
alpha   = atan2(cart(2),cart(1)) %rotation around z-axis
% $$$ if cart(1) < 0
% $$$     alpha = alpha + pi;
% $$$ end
s = sqrt(cart(1)^2+cart(2)^2);
if  0<= cart(1)
    alpha = asin(cart(2)/s)
else
    alpha = pi - asin(cart(2)/s)
end
beta = acos(dot([0 0 1], cart))   %rotation around new x-axis
%theta = acos(cart(3));          %rotation around new x-axis
gamma   = 0; % rotation around new z-axis -- this rotation is not used.

euler = [alpha beta gamma];
basis = euler2basis_vec(euler, [0 1 0]);
%basis
end



function output = euler2basis_vec(euler, vector)

alpha = euler(1);
beta  = euler(2);
gamma = euler(3);

sa = sin(alpha);
ca = cos(alpha);
sb = sin(beta);
cb = cos(beta);
sc = sin(gamma);
cc = cos(gamma);
 
output(1)= sb*sc;
output(2)= cc*sb;
output(3)= cb;

output2(1) = vector(1)*(ca*cc-cb*sa*sc) + vector(2)*(-cb*cc*sa-ca*sc) + vector(3)*sa*sb;
output2(2) = vector(1)*(cc*sa+ca*cb*sc) + vector(2)*(ca*cb*cc-sa*sc)  + vector(3)*(-ca*sb);
output2(3) = vector(1)*sb*sc + vector(2)*cc*sb + vector(3)*cb;

output
output2
% $$$ output(1)= ca*cc-cb*sa*sc;
% $$$ output(2)=-cb*cc*sa-ca*sc;
% $$$ output(3)= sa*sb;
end


% TEST/DeBug
% $$$ 
% $$$ 
% $$$ x = [1;0;0];
% $$$ y = [0;1;0];
% $$$ z = [0;0;1];
% $$$ 
% $$$ R1 = beam2euler_angles(pi/4,0)
% $$$ R2 = beam2euler_angles(0,pi/4)
% $$$ 
% $$$ figure;
% $$$ hold on;
% $$$ plot3([0 x(1)], [0 x(2)], [0 x(3)])
% $$$ plot3([0 y(1)], [0 y(2)], [0 y(3)])
% $$$ plot3([0 z(1)], [0 z(2)], [0 z(3)])
% $$$ plot3(x(1), x(2), x(3), 'b*', 'MarkerSize', 10)
% $$$ plot3(y(1), y(2), y(3), 'b*', 'MarkerSize', 10)
% $$$ plot3(z(1), z(2), z(3), 'b*', 'MarkerSize', 10)
% $$$ 
% $$$ x1 = R1*x;
% $$$ y1 = R1*y;
% $$$ z1 = R1*z;
% $$$ plot3([0 x1(1)], [0 x1(2)], [0 x1(3)], 'r')
% $$$ plot3([0 y1(1)], [0 y1(2)], [0 y1(3)], 'r')
% $$$ plot3([0 z1(1)], [0 z1(2)], [0 z1(3)], 'r')
% $$$ plot3(x1(1), x1(2), x1(3), 'r*', 'MarkerSize', 10)
% $$$ plot3(y1(1), y1(2), y1(3), 'r*', 'MarkerSize', 10)
% $$$ plot3(z1(1), z1(2), z1(3), 'r*', 'MarkerSize', 10)
% $$$ 
% $$$ 
% $$$ x2 = R2*x;
% $$$ y2 = R2*y;
% $$$ z2 = R2*z;
% $$$ plot3([0 x2(1)], [0 x2(2)], [0 x2(3)], 'g')
% $$$ plot3([0 y2(1)], [0 y2(2)], [0 y2(3)], 'g')
% $$$ plot3([0 z2(1)], [0 z2(2)], [0 z2(3)], 'g')
% $$$ plot3(x2(1), x2(2), x2(3), 'g*', 'MarkerSize', 10)
% $$$ plot3(y2(1), y2(2), y2(3), 'g*', 'MarkerSize', 10)
% $$$ plot3(z2(1), z2(2), z2(3), 'g*', 'MarkerSize', 10)
% $$$ 
% $$$ 
% $$$ v1 = R1*z;
% $$$ v2 = R2*v1;
% $$$ plot3([0 v2(1)], [0 v2(2)], [0 v2(3)], 'y')
% $$$ plot3(v2(1), v2(2), v2(3), 'y*', 'MarkerSize', 10)
% $$$ 
% $$$ plot3([0 v1(1)], [0 v1(2)], [0 v1(3)], 'k')
% $$$ plot3(v1(1), v1(2), v1(3), 'k*', 'MarkerSize', 10)
% $$$ 
% $$$ xlabel('X')
% $$$ ylabel('y')
% $$$ axis equal
% $$$ grid on
% $$$ view(3)
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ % $$$ % rotation matrices
% $$$ % $$$ Y = [cos(zx) 0 sin(rot_y); ...
% $$$ % $$$     0 1 0; ...
% $$$ % $$$     -sin(rot_y) 0 cos(rot_y)];
% $$$ % $$$ 
% $$$ % $$$ X = [1 0 0; ...
% $$$ % $$$     0 cos(rot_x) -sin(rot_x); ...
% $$$ % $$$     0 sin(rot_x) cos(rot_x)];
% $$$ % $$$ 
% $$$ % $$$ R = X*Y;  % X and Y commute (X*Y=Y*X)
% $$$ % $$$ 
% $$$ % $$$ euler = R ;%transpose(R*[0;0;1]);
% $$$ % $$$ %return
% $$$ % $$$ % determine euler angles that would give this rotation matrix (www.gregslabaugh.name/publications/euler.pdf)
% $$$ % $$$ if R(3,1) ~= 1 && R(3,1) ~= -1
% $$$ % $$$     % θ
% $$$ % $$$     theta = -asin(R(3,1));
% $$$ % $$$     if theta > pi/2
% $$$ % $$$         theta = pi - theta;
% $$$ % $$$     end
% $$$ % $$$     % ψ
% $$$ % $$$     if cos(theta) > 0
% $$$ % $$$         psi = atan2(R(3,2),  R(3,3));
% $$$ % $$$     else 
% $$$ % $$$         psi = atan2(-R(3,2), -R(3,3));
% $$$ % $$$     end
% $$$ % $$$     % φ
% $$$ % $$$     if cos(theta) > 0
% $$$ % $$$         phi = atan2(R(2,1), R(1,1));
% $$$ % $$$     else
% $$$ % $$$         phi = atan2(-R(2,1), -R(1,1));
% $$$ % $$$     end
% $$$ % $$$ else 
% $$$ % $$$     phi = 0
% $$$ % $$$     if R(3,1) == -1
% $$$ % $$$         theta = pi/2;
% $$$ % $$$         psi   = phi + atan2(R(1,12) , R(1,13));
% $$$ % $$$     else
% $$$ % $$$         theta = -phi/2;
% $$$ % $$$         psi   = -phi + atan2(-R(1,2) , -R(1,3));
% $$$ % $$$     end
% $$$ % $$$ end
% $$$ % $$$ 
% $$$ % $$$ euler = [theta psi phi];
% $$$ 
