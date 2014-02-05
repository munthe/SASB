%  Field II definition of the 32x32 2D matrix transducer
%
%  Parameters:
%   Type:             2D Matrix
%   Center frequency: 3.5 MHz
%   No Elements:      32x32 (35x32) (x y)
%   Bandwidth:        60 %
%   Pitch:            0.3 mm
%   Kerf:             lambda/100 
%
%  Calling: [th, thi] = vermon_35x32 (fs, 'c', 1480, 'name', 'dense_64x64');
%
%  Input:   fs      - Sampling frequency [Hz]
%           c       - Speed of sound     [m/s]
%           plot    - when true, information about the transducer is plotted.
%           name    - Name of transducer type. Possible values are:
%                     vermon_32x32, vermon_32x22, vermon_32x16, vermon_32x6, dense_32x32 and dense_64x64.
%           pitch   - pitch of array. Standard is 300e-6 m.
%
%  Output:  th      - Transducer handle from Field II
%           thi     - Transducer handle info
%
%  Standard values of the input parameters are:
%           fs      - No standard value. Must be set.
%           c       - 1540 m/s
%           name    - vermon_32x32
%

%
%  Version 0.01, 16/09-2010  by MFR -- init version.
%          0.02, 28/09-2010  by MFR -- fs and c are now input vars.
%          0.03, 29/09-2010  by MFR -- enabled is now also an input var.
%          0.04, 02/11-2010  by MFR -- Expanded size to 35x32 by adding three inactive columns.
%          0.05, 10/11-2010  by MFR -- Changed input 'enabled' matrix from 32x32 to 35x32. (No more adding columns)
%          0.06  27/02-2011  by MFR -- Changed impulse response to be a simple sinusoid windoved by a hanning.
%          0.07  07/03-2011  by MFR -- Renamed from 35x32 to 32x32. Removed option to set 'enable'. Added
%                                      output: 'no_elements'.
%          0.08  19/07-2011  by MFR -- Added input: 'conn_plugs'.
%          0.09  06/11-2011  by MFR -- Added Vermon 32x32 impulse response
%          0.10  15/11-2011  by MFR -- Removed input 'conn_plugs', and added intput 'name'.
%          0.11  14/08-2012  by MFR -- Added pitch as input argument.
%

function [th, par] = vermon_32x32 (fs, varargin)

%%  Define the various parameters
if nargin < 1, error('no sampling frequenzy was set.'), end
opt.c  = 1540; %[m/s]
opt.plot = 0;
%opt.conn_plugs = [1 2 3 4 5 6]; % connected plugs
%opt.dense = 0;
%opt.pos = [];
opt.pitch = 300e-6; %[m]
opt.name = 'vermon_32x32';
opt = va_arg(opt, varargin);


%Frequency attributes
par.f0 = 3.5e6;              %  Center frequency            [Hz]
%Internal variable
lambda = opt.c/par.f0;         %  Wave length --Only used internally [m]
par.no_elements = 0;
par.name = opt.name;


% --------------------------------------------------
% Transducer Specific Setup (element positions)
% --------------------------------------------------
if (strcmp(par.name,'vermon_32x32')) % Full SARUS
    par.no_elements_x = 32;
    par.no_elements_y = 32;
    par.no_elements   = par.no_elements_x * par.no_elements_y; 
    no_dead_rows = 3;
    dead_row_idx = [27 18 9];
    
elseif (strcmp(par.name,'vermon_32x22')) % Left Rack + 192 Channels from right side = 704 channels 
    par.no_elements_x = 32;
    par.no_elements_y = 22;
    par.no_elements = par.no_elements_x * par.no_elements_y; 
    no_dead_rows = 2;
    dead_row_idx = [18 9];

elseif (strcmp(par.name,'vermon_32x16')) % Left Rack = 512 Channels
    par.no_elements_x = 32;
    par.no_elements_y = 16;
    par.no_elements = par.no_elements_x * par.no_elements_y; 
    no_dead_rows = 1;
    dead_row_idx = [9];

elseif (strcmp(par.name,'vermon_32x6')) % 1 Plug = 192 Channels
    par.no_elements_x = 32;
    par.no_elements_y = 6;
    par.no_elements = par.no_elements_x * par.no_elements_y; 
    no_dead_rows = 0;
    dead_row_idx = [];
    
    
    % ---------------------------------------------
    % Simulation Probes
    % ---------------------------------------------
elseif (strcmp(par.name,'dense_32x32'))
    par.no_elements_x = 32;
    par.no_elements_y = 32;
    par.no_elements = par.no_elements_x * par.no_elements_y; 
    no_dead_rows = 0;
    dead_row_idx = [];
    
elseif (strcmp(par.name,'dense_64x64'))
    par.no_elements_x = 64;
    par.no_elements_y = 64;
    par.no_elements = par.no_elements_x * par.no_elements_y; 
    no_dead_rows = 0;
    dead_row_idx = [];

else
    error(sprintf([' Unsupported aperture.\n Known apertures are: \n'...
                   '  "vermon_32x6", "vermon_32x16", "vermon_32x22", "vermon_32x32",\n' ...
                   '  "dense_32x32" and "dense_64x64."']));
end


% Enabled Matrix
par.enabled = ones(par.no_elements_x, par.no_elements_y+no_dead_rows);
% Remove dead rows (if any)
par.enabled(:,dead_row_idx) = 0;


%X-dim
par.pitch  = opt.pitch;        %  Element pitch               [m]
par.Nx     = size(par.enabled,1); %  Number of elements in the x-direction
par.kerf_x = lambda/100;    %  Element kerf in x-direction [m]
par.width  = par.pitch-par.kerf_x; % size in x-dimension       [m]
par.no_sub_x = 1;           %  Number of sub-divisions in x-direction of elements
%Y-dim
par.Ny=size(par.enabled,2); %  Number of elements in the y-direction
par.kerf_y = par.kerf_x;    %  Element kerf in y-direction [m]
par.height = par.width;     %  size in y-domension         [m]
par.no_sub_y = 1;           %  Number of sub-divisions in y-direction of elements

%Electronic focus
par.focus = [0 0 0]/1000; %[m]  -- hvad bruges den egentlig til??


%% Create Transducer in Field
%Define the transducer layout
th =  xdc_2d_array (...
    par.Nx, ...
    par.Ny, ...
    par.width, ...
    par.height, ...
    par.kerf_x, ...
    par.kerf_y, ...
    par.enabled, ...
    par.no_sub_x, ...
    par.no_sub_y, ...
    par.focus);


%% Get element positions
info  = xdc_get(th, 'rect');
par.elem_pos = [info(24,1:par.no_sub_x*par.no_sub_y:end)' ...
                info(25,1:par.no_sub_x*par.no_sub_y:end)' ...
                info(26,1:par.no_sub_x*par.no_sub_y:end)']; % [x y z]



%% Pulse-Echo Impulse Response
par.h = impulse_vermon35x32(fs);
par.fs=fs;  % Sampling frequency  [Hz]

[val par.offset] = max(abs(hilbert(par.h)));

% $$$ %use simple excitation pulse instead
% $$$ par.n_pulses = 3;
% $$$ par.h=sin(2*pi*par.f0*(0:1/fs:par.n_pulses/par.f0));
% $$$ par.h = par.h .* hanning(length(par.h))';

%Tell field about the impulse response
xdc_impulse (th, par.h);


