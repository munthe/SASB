function [pos amp] = cfu_cyst_phantom(medium, cyst, silent, seed)
% [pos amp] = cfu_cyst_phantom(medium, cyst, [silent, seed])
% Function that defines a Field 2 cyst phantom.
%
% medium.x    = [x_min x_max] - x-limit of medium in meters
% medium.y    = [y_min y_max] - y-limit of medium in meters
% medium.z    = [z_min z_max] - z-limit of medium in meters
% medium.dens = average number of scatterer per mm^3
%
% cyst.r      = [r_x r_y r_z] - radius in the x, y and z-dimension
% cyst.c      = [c_x c_y c_z] - center coordinate of the cyst 
% cyst.amp    = amplitude of cyst (set to 0 for anechoic, 1 for same amplitude as medium)
%
% Silent:   Set to true to disable printing info to the terminal.
% seed:     Set to an interger between 0 and 2^32-1 to force a seed. Otherwise the seed
%             is random, based on the current time.
%
% By MOFI
% 2013-12-12 v1.0 Init version.
% 2013-12-12 v1.1 Now prints to stdout. Added 'silent' argument.
% 2013-12-12 v1.2 Now uses a local random seed stream. Added 'seed' argument.
%


if nargin < 3, silent = false; end
if nargin < 4, seed   = []; end

% TODO: Add argument verification



% From input structs to vars
x_min = medium.x(1);
x_max = medium.x(2);
y_min = medium.y(1);
y_max = medium.y(2);
z_min = medium.z(1);
z_max = medium.z(2);

r_x    = cyst.r(1);
r_y    = cyst.r(2);
r_z    = cyst.r(3);

c_x    = cyst.c(1);
c_y    = cyst.c(2);
c_z    = cyst.c(3);

% calc helper vars
mm_to_m = 1e-3;
m_to_mm = 1e3;

if y_min == y_max
    area = abs((x_max-x_min) * (z_max-z_min));
    N    = ceil(area*(m_to_mm)^2 *medium.dens);
    if ~silent, fprintf('Cyst phantom: medium area is %.1fmm^2.\n', area*(m_to_mm)^2);end
else
    volume  = abs((x_max-x_min) * (y_max-y_min) * (z_max-z_min));
    N   = ceil(volume*(m_to_mm)^3  * medium.dens);
    if ~silent, fprintf('Cyst phantom: medium volume is %.1fmm^3.\n', volume*(m_to_mm)^3);end
end

% random seed 
seed   = sum(clock);
seed   = mod(seed*1e7, 2^32);
stream = RandStream('mt19937ar','Seed',seed);

% make medium coordinates
pos_x = rand(stream,N,1)*(x_max-x_min) + x_min;
pos_y = rand(stream,N,1)*(y_max-y_min) + y_min;
pos_z = rand(stream,N,1)*(z_max-z_min) + z_min;
% make medium random amplitudes
amp   = randn(stream,N,1);


% Diskriminate between inside and outside of cyst 
d = (((pos_x - c_x)./r_x).^2 + ((pos_y - c_y)./r_y).^2 +((pos_z - c_z)./r_z).^2).^(1/2);
d = (abs(d)<1);

% Correct the amplitude of the scatterer inside the cyst (or remove them if amp==0).
if cyst.amp == 0
    pos_x(d) = [];
    pos_y(d) = [];
    pos_z(d) = [];
    amp(d)   = [];
else
    amp(d)   = amp(d).*cyst.amp;
end

if ~silent, fprintf('Cyst phantom: Number of scatterer defined %i.\n', length(amp));end

% output position matrix
pos = [pos_x pos_y pos_z];
