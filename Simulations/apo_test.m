% Simulates imaging with different numbers of indexed values of Hanning
% apodization vector.
%
% Assumes that field_init has been run
%

%% Initilize
initilize

%% Apodization

% Make the apodization vector
% apo = ones(1,N_active);
% apo = hanning(N_active)';

trans_focus = input('Transmit focus [mm]: ')/1000;
receive_focus = input('Receive focus [mm]: ')/1000;
apo_steps = [1,2,4,8,16];
num_plots = length(apo_steps);
figure
for j = 1:num_plots
    subplot(1,num_plots,j)
    hold on
    
    % Calculate indexed apodization vector
    apo = ceil(hanning(N_active)'.*apo_steps(j))./apo_steps(j);

    sesr
    mk_img
    
    title(strcat(num2str(apo_steps(j)),' apodization levels'))
    if j==1;
        axis on
        ylabel('Axial distance [mm]') 
        xlabel('Lateral distance [mm]')
    end
    hold off
end