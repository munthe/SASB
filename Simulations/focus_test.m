% Simulates imaging with different foci and apodization vectors.
%
% Assumes that field_init has been run
%

%% Initilize
initilize

%% Apodization

% Make the apodization vector
% apo = ones(1,N_active);
apo = hanning(N_active)';

[trans_foci] = input('Transmit foci vector: ');
[receive_foci] = input('Receive foci vector (should be of same length as transmit): ');
num_plots = length(trans_foci);
figure
for j = 1:num_plots
    subplot(1,num_plots,j)
    hold on

    % Set focus
    trans_focus = trans_foci(j)/1000;
    receive_focus = receive_foci(j)/1000;

    sesr
    mk_img
    
    scatter([-10,10],[trans_foci(j),receive_foci(j)],100,'Xr')
    if j==1;
        axis on
        ylabel('Axial distance [mm]') 
        xlabel('Lateral distance [mm]')
    end
    hold off
end