% $Id: psf_sasb_8804.m,v 1.15 2011-07-17 21:47:07 jmh Exp $
if isunix
  addpath('/usr/local/home/munk/programming/matlab/toolbox/Field_II_linux_7_10');
  addpath('../src');
  folder = [tempdir 'orbit-munk/'];
else
  addpath('z:/programming/matlab/aperture/src');
  folder = 'c:/tmp/';
end

% For cluster usage
silent = 0;

if ~exist('use_field','var')
  use_field = 1;
end

% Interpolation
interpolation = 'spline';

% Decimation
downsmple = 1;

% Transducer - 8804
f0=7e6;    % Transducer center frequency [Hz]
fs=120e6;  % Sampling frequency
c=1540;    % Speed of sound [m/s]

n_elements = 192;
pitch = 208e-6;
array_width = 40e-3;

kerf = 0.035 / 1000;
height=4.5/1000;
width = pitch - kerf;
r_focus = 20.0/1000; % Elevation focus
lambda = c/f0;

e_focus_depth = 10/1000;

% Create an image
n_lines = 192;
dx = pitch/2;

% Excitation and impulse response
n_impulse_cycles = 3;
n_excitation_cycles = 2;

n_active_elements = 64; % Used for calculation resolution cell
n_sub_x = 10;
n_sub_y = 1;
n_sub_y_focused = 10;

% Anonymous function for logarithmic compression
log_compression = @(x) 20*log10(x/max(x(:)));

if use_field
  try
    field_end
  catch
    disp('Error in resetting FieldII\n')
  end

  field_init(-1);

  set_field('c',c);   % Set speed of sound
  set_field('fs',fs); % Set the sampling frequency

  % Electronic focus
  focus=[0 0 e_focus_depth];    % Changed for each line

  % Generate aperture for emission
  xmt_aperture = xdc_focused_array(n_elements, width, height, kerf,...
                                   r_focus, n_sub_x,...
                                   n_sub_y_focused, focus);

  % Set excitation of the emit aperture
  excitation=sin(2*pi*f0*(0:1/fs:n_excitation_cycles/f0));
  xdc_excitation(xmt_aperture, excitation);

  % Set the impulse response
  impulse_response = sin(2*pi*f0*(0:1/fs:n_impulse_cycles/f0));
  impulse_response = impulse_response.* hanning(max(size(impulse_response)))';
  xdc_impulse(xmt_aperture, impulse_response);

  % Generate aperture for reception
  rcv_aperture = xdc_focused_array (n_elements, width, height, ...
                                    kerf, r_focus, n_sub_x, ...
                                    n_sub_y_focused, focus);

  % Set the impulse response for the receive aperture
  xdc_impulse (rcv_aperture, impulse_response);

end
% Phantom
phantom_positions = (10:5:60)./1000;
phantom_positions = [zeros(size(phantom_positions)) ; ...
                    zeros(size(phantom_positions)) ; phantom_positions]';
phantom_amplitudes = ones(size(phantom_positions,1),1).*1.0;

tmp = 0.01+max(sqrt(sum(phantom_positions.^2,2))) - ...
      min(sqrt(sum(phantom_positions.^2,2))); % 

n_rf_samples = ceil(2* tmp/c *fs+1);

% Aperture position
xpos = -((n_elements * pitch) - kerf)/2 + 0.5*width : pitch ...
       : ((n_elements * pitch) - kerf)/2 - 0.5*width;
ypos = zeros(1,n_elements);
zpos = zeros(1,n_elements);
pos = [xpos ; ypos ; zpos]';

% Line positions
offset_x = 0;
wx = (n_lines-1)/2 + offset_x;
line_x = ([0:n_lines-1] - wx) * dx;

origin =   [line_x ; zeros(1,n_lines) ...
                   ; 0/1000.*ones(1,n_lines)]';

w_hamming = @(n)  0.46*cos(pi*n) + 0.54;
w_apodization = w_hamming;

for line_no=1:n_lines
  disp(sprintf('Field II: Line %d',line_no))
  if exist(sprintf('%sline%d.mat',folder,line_no),'file')
    continue
  else
    system(sprintf('touch %sline%d.mat',folder,line_no));
  end
  
  xdc_focus_times(xmt_aperture, 0, zeros(1, n_elements));
  xdc_focus_times(rcv_aperture, 0, zeros(1, n_elements));

  % Focus
  focus = [origin(line_no,1) 0 e_focus_depth];

  xdc_center_focus(xmt_aperture, [focus(1) 0 0]); 
  
  xdc_focus(xmt_aperture,0,focus);

  % Apodization
  apodization = zeros(1,n_elements);
  index_norm = abs(xpos - origin(line_no,1)) / (n_active_elements*pitch/2);
  apodization(index_norm < 1) = w_apodization(index_norm(index_norm < 1));

  xdc_apodization(xmt_aperture, 0, apodization);
  xdc_apodization(rcv_aperture, 0, ones(1,n_elements));
  
  [rf_data, tstart] = calc_scat_multi(xmt_aperture, rcv_aperture, ...
                                      phantom_positions, phantom_amplitudes);

  rf_data = rf_data(:,find(apodization>0));
  
  
  % Correct for length of responses
  [val index] = max(abs(hilbert(conv(conv(impulse_response,excitation), ...
                                     impulse_response))));
  tstart = tstart - (index-1)/(fs);
  
  % Correct for elevation focus
  v = atan2((height/2), r_focus);
  delay = ((r_focus / cos(v)) - r_focus)/c*4;
  tstart = tstart - delay;
  
  % Store RF data
  save(sprintf('%sline%d.mat',folder,line_no),'rf_data');
  save(sprintf('%sline%d.mat',folder,line_no),'tstart','-APPEND');
  
  % Beam formation with Field II
  xdc_center_focus(rcv_aperture,[focus(1) 0 0]);
  xdc_dynamic_focus(rcv_aperture,-100,0,0);
  xdc_apodization(rcv_aperture, 0, apodization);

  [rf_data, tstart] = calc_scat(xmt_aperture, rcv_aperture,...
                                phantom_positions, phantom_amplitudes);

  % Correct for length of responses
  [val index] = max(abs(hilbert(conv(conv(impulse_response,excitation), ...
                                     impulse_response))));
  tstart = tstart - (index-1)/(fs);
  
  % Correct for elevation focus
  v = atan2((height/2), r_focus);
  delay = ((r_focus / cos(v)) - r_focus)/c*4;
  tstart = tstart - delay;
  
  field_times(line_no) = tstart;
  field_data_store(1:size(rf_data,1),line_no) = rf_data(:,:); 
end

if exist(sprintf('%sfield_times.mat',folder),'file')
  pp = open(sprintf('%sfield_times.mat',folder));
  field_times = pp.field_times;
  min_sample = pp.min_sample;
  n_rf_samples = pp.n_rf_samples;
else
  min_sample = min(field_times)*fs;
  max_sample = max(field_times)*fs;
  [n,m]=size(field_data_store(:,:));
  n_rf_samples = n + (max_sample-min_sample);
  
  env_field = zeros(round(n_rf_samples),n_lines);
  
  for line_no=1:n_lines
    rf_env = abs(hilbert([zeros(round(field_times(line_no)*fs-min_sample),1); ...
                        field_data_store(:,line_no)]));
    env_field(1:max(size(rf_env)),line_no)=rf_env;
  end
  
  env_field = log_compression(env_field); % Avoid NaN
  
  if ~silent
    figure(1);
    xlim = [line_x(1) line_x(end)]*1000;
    ylim = 1000*min_sample*c/(2*fs) + 1000.*[0 n_rf_samples]*c/(2*fs);
    imagesc(xlim,ylim, env_field,[-60 0]);
    axis('image')
    colormap(gray)
    xlabel('Vertical position x [mm]')
    ylabel('Depth z [mm]')
    title('Field II')
  end

  save(sprintf('%sfield_times.mat',folder),'field_times');
  save(sprintf('%sfield_times.mat',folder),'n_rf_samples','-APPEND');
  save(sprintf('%sfield_times.mat',folder),'min_sample','-APPEND');
end

if silent
  exit;
end

if downsmple
  n_rf_samples = round(n_rf_samples / 3);
  fs = fs / 3;
end

% Origin of lines
origin(:,3) =  min_sample*c/(2*fs).*ones(n_lines,1);

% Lines
direction = repmat([0 0 1],[n_lines 1]);
dr = repmat(c/(2*fs),[n_lines 1]);

% Length
line_length = c*((n_rf_samples-min_sample)/(2*fs)).*ones(n_lines,1);

bft2_n_rf_samples = floor(n_rf_samples);

tic

globals = bft3_system('fs',fs,'c',c);

xmt_aperture = bft3_aperture('type','linear_array',...
                             'pitch',pitch, 'n_elements',n_elements);

rcv_aperture = bft3_aperture('type','linear_array',...
                             'pitch',kerf, 'n_elements',n_active_elements);

bft_lines = [];
xmt_apodizations = [];
rcv_apodizations = [];

apo = ones(size(xmt_aperture.pos,1),1);
apotimes = [0] / 1000;

for line_no=1:n_lines
  % XMT apodization
  aporef = [origin(line_no,1:2) 0.0];
  tmp =  bft3_apodization(xmt_aperture, aporef, apotimes, apo);
  xmt_apodizations = [xmt_apodizations tmp];
end

apodization = hamming(n_active_elements);
for line_no=1:n_lines
  % RCV apodization
  tmp = bft3_apodization(rcv_aperture, aporef, apotimes, apodization);
  rcv_apodizations = [rcv_apodizations tmp];
  
  % Line
  bft_lines = [bft_lines ...
               bft3_line(origin(line_no,:), direction(line_no,:),...
                         dr(line_no), line_length(line_no))];
end

toc

if ~exist('bf_image_bft3')
  bf_image_bft3 = zeros(size(bft_lines(1).pos,1), n_lines);

  rcv_pos = rcv_aperture.pos;
  tic
    for line_no=1:n_lines
      disp(['BFT III: emission no: ' num2str(line_no)]);

      % Set positions to active elements
      index_norm = abs(xpos - origin(line_no,1)) / (n_active_elements*pitch/2);
      rcv_aperture.pos = xmt_aperture.pos(find(index_norm < 1),:);
      
      bft_image = bft3_image(xmt_aperture, rcv_aperture, ...
                             xmt_apodizations(line_no),...
                             rcv_apodizations(line_no),...
                             bft_lines(line_no));

      bft_image.interp = interpolation;
      bft_image.nthreads = int32(1);
      
      focus = [origin(line_no,1) 0  e_focus_depth];
      
      xmt_aperture.center_focus = [origin(line_no,1) 0 0]; 

      % Virtual source in xmt and rcv
      xmt_aperture.focus = focus;
      rcv_aperture.focus = focus;
      
      % Read file
      pp = open(sprintf('%sline%d.mat',folder,line_no));
      if downsmple
        rf_data = downsample(pp.rf_data,3);
      else
        rf_data = pp.rf_data;
      end
      
      tstart = pp.tstart;
      
      bf_image_bft3(:,line_no) = ...
          bft_image.beamform(rf_data, tstart, uint32(line_no));
    end
  toc
end

env = abs(hilbert(bf_image_bft3(:,:)));
env = log_compression(env);

figure(2);

xlim = [line_x(1) line_x(end)]*1000;
ylim = [origin(1,3) origin(1,3)+line_length(1)]*1000;
imagesc(xlim, ylim, env,[-60 0]);

colormap(gray)
xlabel('Vertical position x [mm]')
ylabel('Depth z [mm]')
title('BFT3')
axis('image')
set(gca,'YLim',ylim)

clear bft_image xmt_apodizations rcv_apodization;

disp('BFT3 setup')
tic

sasb_xmt_aperture = bft3_aperture('pos',[line_x ; line_x*0 ; line_x*0]');
sasb_rcv_aperture = bft3_aperture('pos',[line_x(1) 0 0]);

xmt_apodizations = [];
rcv_apodizations = [];

apo = hamming(size(sasb_xmt_aperture.pos,1));
apotimes = [0] / 1000;

for i=1:n_lines
  % XMT apodization
  aporef = [origin(i,1:2) e_focus_depth];
  tmp =  bft3_apodization(sasb_xmt_aperture, aporef, apotimes, apo);
  tmp.parametric  = true;
  tmp.dynamic = true;
  tmp.f = 2;
  xmt_apodizations = [xmt_apodizations tmp];
end

apo = ones(1,1);
for i=1:n_lines
  % RCV apodization
  aporef = [origin(i,1:2) e_focus_depth];
  tmp = bft3_apodization(sasb_rcv_aperture, aporef, apotimes, apo);
  tmp.parametric  = false;
  tmp.dynamic = true;
  rcv_apodizations = [rcv_apodizations tmp];
end

my_sasb_image = zeros(size(bf_image_bft3));

bft_image = bft3_image(sasb_xmt_aperture, sasb_rcv_aperture, ...
                       xmt_apodizations, rcv_apodizations,...
                       bft_lines);
bft_image.interp = interpolation;
bft_image.nthreads = int32(4);
toc

disp('BFT3 beam formation')
tic

line_indices = 1:n_lines;
  
for line_index=1:length(line_indices)
  i = line_indices(line_index);
%  disp(['BFT III: SASB line no: ' num2str(i)]);
  sasb_xmt_aperture.center_focus = [line_x(i) 0 0];
  sasb_xmt_aperture.focus = [line_x(i) 0 e_focus_depth];

  sasb_rcv_aperture.pos = [line_x(i) 0 0];
  sasb_rcv_aperture.focus = [line_x(i) 0 e_focus_depth];

  my_image = bft_image.beamform(bf_image_bft3(:,i),...
                                min_sample/fs,uint32(i));

  my_sasb_image = my_sasb_image + my_image;
end
toc

env_sasb = abs(hilbert(my_sasb_image));
env_sasb = log_compression(env_sasb);

figure(3);
imagesc(xlim, [origin(1,3) origin(1,3)+line_length(1) ...
                   ]*1000,env_sasb,[-60 0]);
colormap(gray)
xlabel('Vertical position x [mm]')
ylabel('Depth z [mm]')
title('BFT3 - SASB')
axis('image')
set(gca,'YLim',ylim)

