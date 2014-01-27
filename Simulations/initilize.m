% Initilized values and transducers
%

%%  Initilize values

f0=3e6;                 % Transducer center frequency [Hz]
fs=100e6;               % Sampling frequency [Hz]
c=1540;                 % Speed of sound [m/s]
lambda=c/f0;            % Wavelength [m]
width=lambda;           % Width of element
element_height=5/1000;  % Height of element [m]
kerf=0.1/1000;          % Kerf [m]
focus=[0 0 00]/1000;    % Initial focal point [m]
N_elements=128;         % Number of physical elements
N_active=64;            % Number of active elements
% xmit_N_active=128;      % Number of active transmit elements for constant F#
% rec_N_active=128;	    % Number of active receive elements for constant F#

% Set the sampling frequency 
set_sampling(fs);

%   Load the computer phantom
[phantom_positions, phantom_amplitudes] = pts_pha;

%   Do linear array imaging
no_lines=20;                         %  Number of lines in image
image_width=20/1000;                 %  Size of image sector
d_x=image_width/no_lines;            %  Increment for image

%% Setup transducer

%  Generate aperture for emission
emit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 1,focus);

%  Set the impulse response and excitation of the emit aperture
impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (emit_aperture, impulse_response);

excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (emit_aperture, excitation);   
    
%  Generate aperture for reception
receive_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 1, 1,focus);

%  Set the impulse response for the receive aperture
xdc_impulse (receive_aperture, impulse_response);
