% cfu_plot_field_xdc.m
% This script plots the transducer layout and indicates each element's current apodization.
% It is inspired by a script from the Field2 manual, but this version actually functions.
%
%  xdc_plot (th)
%
% By Morten F. Rasmussen, 
% Version 1.0  2011-03-10 Init version
% Version 1.1  2011-11-09 Remove edge color
%

function cfu_plot_field_xdc (th)

data            = xdc_get(th,'rect');
if isempty(data)
    data            = xdc_get(th,'lin');
    if isempty(data)
        data            = xdc_get(th,'tri');
    end
end


apo             = xdc_get(th, 'apo');
apo             = apo(2:end); % remove time
no_mat_elements = size(data,2);
no_elements     = size(apo,1);
no_mat_per_elem = no_mat_elements/no_elements;
% expand element apo to math-element apo
mat_apo         = repmat(apo,1,no_mat_per_elem)';
mat_apo         = mat_apo(:);

if no_mat_elements > 10000,
    fprintf('Warning: this will loop through %i mathematical elements.\n',no_mat_elements);
    disp('Are you sure you want to continue?');
    disp('(Hit enter to continue or ctr-c to quit)');
    pause;
    disp('Continuing...');
end


hold on;
% Plot elements
elm_idx=1;
for elem_i=1:no_mat_elements,
    if(rem(elem_i,1000) == 1)
        drawnow
    end
    x = [data(11,elem_i), data(20,elem_i);  data(14,elem_i), data(17,elem_i)]*1000;
    y = [data(12,elem_i), data(21,elem_i);  data(15,elem_i), data(18,elem_i)]*1000;
    z = [data(13,elem_i), data(22,elem_i);  data(16,elem_i), data(19,elem_i)]*1000;
    %fprintf('elem_i: %i  apo: %4.2f\n', elem_i, apo(elem_i));
    c = mat_apo(elem_i)*ones(2,2);
    %c = apo(elm_idx)*ones(2,2);
    %surf(x,y,z,c, 'EdgeColor', 'None');
    surf(x,y,z,c);
    elm_idx= elm_idx+1;
end
cbh = colorbar;
%view(3)
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')
grid on;
axis('image');
hold off;
% $$$ set(gca, 'XDir', 'reverse')
% $$$ set(gca, 'YDir', 'reverse')

