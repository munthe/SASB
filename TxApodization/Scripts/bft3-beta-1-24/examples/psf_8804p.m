%> @file psf_8804p.m
%> @brief Beamformation example
%>
%> It is assumed that FieldII is installed and the location
%> included in the path. If not correct the addpath below.
%>
% $Id: psf_8804p.m,v 1.14 2011-07-17 14:32:08 jmh Exp $
if isunix
  [dummy user] = system('echo $USER | tr -d "\n"');
%  addpath(sprintf(['/home/%s/programming/matlab/toolbox/' ...
%                   'Field_II_linux_7_10/'],user));
  addpath(sprintf('/usr/local/home/%s/programming/matlab/toolbox/Field_II_linux_7_10/',user));
  folder = [tempdir 'orbit-' user '/'];
  system(sprintf('mkdir -p %s',folder));
else
  error('Your operating system came from a monopoly: Upgrade to *nix');
end

addpath(['../src']);

interp = 'linear';

% Transducer - 8804
f0=7e6;   % Transducer center frequency [Hz], lower to 5 Mhz for aoi > 15
fs=120e6; % Sampling frequency
c=1540;   % Speed of sound [m/s]

n_elements = 192;
pitch = 208e-6;
array_width = 40e-3;

kerf = 0.035 / 1000;
height=4.5/1000;
width = pitch - kerf;
r_focus = 20.0/1000;   % Elevation focus
lambda = c/f0;

n_active_elements = 64;
n_sub_x = 1;
n_sub_y = 10;
n_sub_y_focused = 10;

e_focus =  40/1000;    % Electronic focus
depth   = 30/1000;

fixed_rcv_focus = 1;

focus   = [0 0 e_focus];    % Changed for each line

% Anonymous function for logarithmic compression, use imagesc to
% truncate to [-dBmax 0] dB
log_compression = @(x) 20*log10(x/max(x(:)));
dBmax = 60;

% Create an image
n_lines=50;         % Number of A-lines in image
dx=pitch;           % Increment for image

try
  field_end
catch
  disp('Error in resetting FieldII')
end

field_init(-1);

set_field('c',c);   % Set speed of sound
set_field('fs',fs); % Set the sampling frequency

% Generate aperture for emission
xmt_aperture = xdc_linear_array(n_elements, width, height, kerf,...
                                n_sub_x, n_sub_y,focus);

n_excitation_cycles = 1;
n_impulse_cycles = 2;

% Set excitation of the emit aperture
excitation=sin(2*pi*f0*(0:1/fs:n_excitation_cycles/f0));
xdc_excitation (xmt_aperture, excitation);

% Set the impulse response
impulse_response = sin(2*pi*f0*(0:1/fs:n_impulse_cycles/f0));
impulse_response = impulse_response.* hanning(max(size(impulse_response)))';


xdc_impulse (xmt_aperture, impulse_response);

% Generate aperture for reception
rcv_aperture = xdc_linear_array (n_elements, width, height, ...
                                 kerf, n_sub_x, n_sub_y, focus);

% Set the impulse response for the receive aperture
xdc_impulse (rcv_aperture, impulse_response);

% Phantom
phantom_positions = [0 0 depth];
phantom_amplitudes = ones(size(phantom_positions,1),1).*1.0;

% Aperture position
xpos = -((n_elements * pitch) - kerf)/2 + 0.5*width : pitch ...
       : ((n_elements * pitch) - kerf)/2 - 0.5*width;

% Line positions
wx = (n_lines-1)/2;
line_x = ([0:n_lines-1] - wx) * dx;

origin =   [line_x ; zeros(1,n_lines) ...
                   ; 0/1000.*ones(1,n_lines)]';

w_hamming = @(n)  0.46*cos(pi*n) + 0.54;
w_hanning = @(n)  0.50*cos(pi*n) + 0.50;
w_boxcar  = @(n)  1.0;

w_apodization = w_hamming;

disp(sprintf(['Field II: Computing Point Spread Function (PSP) at %d ' ...
              'mm'],1000*depth))

for line_no=1:n_lines
  if exist(sprintf('%sline%d.mat',folder,line_no),'file')
    continue;
  else
    dummy = 0;
    save(sprintf('%sline%d.mat',folder,line_no),'dummy');
  end
  disp(sprintf('Field II: line %d',line_no));
  
  xdc_focus_times(xmt_aperture, 0, zeros(1, n_elements));
  xdc_focus_times(rcv_aperture, 0, zeros(1, n_elements));

  apodization = zeros(1,n_elements);
  index_norm = abs(xpos - origin(line_no,1)) / (n_active_elements*pitch/2);
  apodization(index_norm < 1) = w_apodization(index_norm(index_norm < 1));
  
  % Focus
  focus = [origin(line_no,1) 0 e_focus];

  xdc_center_focus(xmt_aperture, [focus(1) 0 0]);
  
  xdc_focus(xmt_aperture,0,focus);

  xdc_apodization(xmt_aperture, 0, apodization);

  xdc_apodization(rcv_aperture, 0, ones(1,n_elements));
  

  [rf_data, tstart] = calc_scat_multi(xmt_aperture, rcv_aperture, ...
                                      phantom_positions, phantom_amplitudes);

  % Correct for length of responses
  [val index] = max(abs(hilbert(conv(conv(impulse_response,excitation), ...
                                     impulse_response))));
  tstart = tstart - (index-1)/(fs);
  
  rf_times(line_no) = tstart;
  
  % Store RF data
  save(sprintf('%sline%d.mat',folder,line_no),'tstart');
  save(sprintf('%sline%d.mat',folder,line_no),'rf_data','-APPEND');

  rf_data_size(line_no) = size(rf_data,1);
  rf_data_store(1:size(rf_data,1),:,line_no) = rf_data(:,:);
  
  % Beamf formation with Field II - dynamic focus
  if (~fixed_rcv_focus)
    xdc_center_focus(rcv_aperture,[focus(1) 0 0]);
    xdc_dynamic_focus(rcv_aperture,-100,0,0);
  else
    % Beam formation with Field II - fixed receive focus
    xdc_focus(rcv_aperture, 0,focus);
  end
  
  xdc_apodization(rcv_aperture, 0, apodization);

  [rf_data, tstart] = calc_scat(xmt_aperture, rcv_aperture,...
                                phantom_positions, phantom_amplitudes);

  % Correct for length of responses
  [val index] = max(abs(hilbert(conv(conv(impulse_response,excitation), ...
                                     impulse_response))));
  tstart = tstart - (index-1)/(fs);
  
  field_times(line_no) = tstart;
  field_data_store(1:size(rf_data,1),line_no) = rf_data(:,:);

  field_tstart = tstart;
  field_rf_data = rf_data;
  save(sprintf('%sline%d.mat',folder,line_no),'field_tstart','-APPEND');
  save(sprintf('%sline%d.mat',folder,line_no),'field_rf_data','-APPEND');
end

if ~exist('field_data_store','var')
  for line_no=1:n_lines
    fstruct = open(sprintf('%sline%d.mat',folder,line_no));

    %Field II
    field_data_store(1:size(fstruct.field_rf_data,1),line_no) = ...
        fstruct.field_rf_data(:,:);
    field_times(line_no) = fstruct.field_tstart;

    % BFT
    rf_data_store(1:size(fstruct.rf_data,1),:,line_no) = fstruct.rf_data(:,:);
    rf_times(line_no) = fstruct.tstart;
    rf_data_size(line_no) = size(fstruct.rf_data,1);
  end
end

min_sample = min(field_times)*fs;
max_sample = max(field_times)*fs;
[n,m]=size(field_data_store(:,:));
field_n_rf_samples = ceil(n + (max_sample-min_sample));

n_rf_samples = size(rf_data_store,1);

field_data = zeros(field_n_rf_samples,n_lines);

for i=1:n_lines
  npad = round(field_times(i)*fs-min_sample);
  field_data(1:(size(field_data_store,1)+npad),i) = [zeros(npad,1); ...
                      field_data_store(:,i)];
end

env_field = abs(hilbert(field_data));
env_field = log_compression(env_field);

fh1 = figure(1);
subplot(3,1,1);
xlims = [-n_lines*pitch/2 n_lines*pitch/2]*1000;
ylims = 1000*min_sample*c/(2*fs) + 1000.*[0 field_n_rf_samples]*c/(2*fs);
imagesc(xlims,ylims, env_field,[-dBmax 0]);
axis('image')
colormap(gray)
xlabel('Lateral distance x [mm]')
ylabel('Depth z [mm]')
title('Field II');

fh2 = figure(2);
subplot(3,1,1);
xs = 1000*meshgrid(line_x, env_field(:,1));
ys = 1000*linspace(min_sample*c/(2*fs),...
                   (min_sample+field_n_rf_samples)*c/(2*fs),field_n_rf_samples);
ys = ndgrid(ys,env_field(1,:));

dBstep = -6;
[tmp, ah] = contour(xs,ys,env_field,dBstep.*[0:1:10]);
xlabel('Lateral distance x [mm]')
ylabel('Time [mm]')
ytick = -6.*[10:-1:0];
yticklabel = num2str(ytick');
yticklabel = [yticklabel, repmat(str2mat(' dB'), [size(yticklabel,1) 1])];
hbar = colorbar('YTick',ytick,'YTickLabel',yticklabel);
axis('ij')
title('Field II dynamic receive focusing');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BFT III
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lines
direction = [0 0 1];
dr = c/(2*fs);
line_length = c*(field_n_rf_samples/(2*fs));

% Change origin of lines to start at first non-zero pixel
origin(:,3) =  min_sample*c/(2*fs).*ones(n_lines,1);

globals = bft3_system('fs',fs,'c',c);

rcv_aperture = bft3_aperture('type','linear_array',...
                             'pitch',pitch, 'n_elements',n_elements);

xmt_aperture = bft3_aperture('type','linear_array',...
                             'pitch',pitch, 'n_elements',n_elements);

xmt_apodizations = [];
rcv_apodizations = [];

apo = ones(n_elements,1);
apotimes = [0] / 1000;
aporef = [0 0 0];

for i=1:n_lines
  % XMT apodization
  tmp =  bft3_apodization(xmt_aperture, aporef, apotimes, apo);
  tmp.parametric = false;
  xmt_apodizations = [xmt_apodizations tmp];
end

for line_no=1:n_lines
  % RCV apodization
  apodization = zeros(1,n_elements);
  index_norm = abs(xpos - origin(line_no,1)) / (n_active_elements*pitch/2);
  apodization(index_norm < 1) = w_apodization(index_norm(index_norm < 1));

  tmp = bft3_apodization(rcv_aperture, aporef, apotimes, apodization');
  tmp.parametric = true;
  rcv_apodizations = [rcv_apodizations tmp];
  
end

% Create multiple lines in one go
bft_lines = bft3_lines(origin,direction,dr,line_length);

bft_n_rf_samples = size(bft_lines(1).pos,1);

bf_image_bft3 = zeros(bft_n_rf_samples, n_lines);

disp('Fixed XMT focus')
disp(sprintf('Computing Point Spread Function (PSP) at %d mm',1000*depth))
for i=1:n_lines
  bft_image = bft3_image(xmt_aperture, rcv_aperture, ...
                         xmt_apodizations(i), rcv_apodizations(i),...
                         bft_lines(i));
  bft_image.interp = interp;

  focus = [origin(i,1) 0 e_focus];

  xmt_aperture.center_focus = [focus(1) 0 0]; 

  % Transmit focus
  xmt_aperture.focus = focus;
  
  % Dynamic receive focusing, so empty receive focus
  rcv_aperture.focus = [];

  my_image = bft_image.beamform(rf_data_store(1:rf_data_size(i),:,i),...
                                rf_times(i),uint32(i));
  bf_image_bft3(:,i) = my_image(:);
end

env = abs(hilbert(bf_image_bft3(:,:)));
env = log_compression(env);

figure(fh1)
subplot(3,1,2);
imagesc(xlims, [bft_lines(1).pos(1,3) bft_lines(1).pos(end,3)]*1000,...
        env,[-dBmax 0]);
colormap(gray)
xlabel('Lateral distance x [mm]')
ylabel('Depth z [mm]')
axis('image')
set(gca,'YLim',ylim)
title('BFTIII - dynamic receive focusing');

figure(fh2)
subplot(3,1,2);
xs = 1000*meshgrid(line_x, env(:,1));
ys = 1000*bft_lines(1).pos(:,3);
ys = ndgrid(ys,env(1,:));

contour(xs,ys,env,dBstep.*[0:1:10]);
xlabel('Lateral distance x [mm]')
ylabel('Time [mm]')
hbar = colorbar('YTick',ytick,'YTickLabel',yticklabel);
axis('ij')
title('BFTIII - dynamic receive focusing');

disp(sprintf('Computing Point Spread Function (PSP) at %d mm',1000*depth))
for i=1:n_lines
  bft_image = bft3_image(xmt_aperture, rcv_aperture, ...
                         xmt_apodizations(i), rcv_apodizations(i),...
                         bft_lines(i));
  bft_image.interp = interp;

  focus = [origin(i,1) 0 e_focus];

  xmt_aperture.center_focus = [focus(1) 0 0]; 

  xmt_aperture.focus = focus;
  rcv_aperture.focus = focus;

  my_image = bft_image.beamform(rf_data_store(1:rf_data_size(i),:,i),...
                                rf_times(i),uint32(i));
  bf_image_bft3(:,i) = my_image(:);
end

env = abs(hilbert(bf_image_bft3(:,:)));
env = log_compression(env);

figure(fh1)
subplot(3,1,3);
imagesc(xlims, [bft_lines(1).pos(1,3) bft_lines(1).pos(end,3)]*1000,...
        env,[-dBmax 0]);
colormap(gray)
xlabel('Lateral distance x [mm]')
ylabel('Depth z [mm]')
axis('image')
set(gca,'YLim',ylim)
title('BFTIII - fixed receive focusing');

figure(fh2)
subplot(3,1,3);
contour(xs,ys,env,dBstep.*[0:1:10]);
xlabel('Lateral distance x [mm]')
ylabel('Time [mm]')
hbar = colorbar('YTick',ytick,'YTickLabel',yticklabel);
axis('ij')
title('BFTIII - fixed receive focusing');
