%
% win = cfu_window (win_name, n, a)
%
% win_name: 'Hanning', 'Hamming', 'Blackman', 'Tukey'
% n: [0 1]. Sample points of window output. Normalised distance to center of window.
%           0 is center and 1 is edge.
% a: parameter to window function (Only for Blackman and Tukey),
%
% 2011-09-19 Init Version. MFR
% 2012-08-28 Added sanity check on input arguments. MFR
% 2012-08-29 Fixed bug in Blackman. MFR
% 2012-11-20 Allow n<0 and n>1 (force to 0 respectively 1).
% 2013-09-16 Renamed to from mfr_window to cfu_window. MFR
%


function win = cfu_window (win_name, n, a)
%% Sanity check
if isempty(n),                         error('The sample point vector, n, is empty.');end
if sum(isnan(n)),                      error('The sample point vector, n, contains NaN');end
if (max(n(:)) > 1) || (min(n(:)) < 0)
    n(n>1)=1;
    n(n<0)=0;
end 

if (nargin >2 && isnan(a)),            error('a contains NaN');end
if ~ischar(win_name),                  error('Input argument "win_name" must be a string.');end


%% Actual work
switch win_name
  case 'Hanning'
    win = 0.50*cos(pi*n) + 0.50;
  case 'Hamming' 
    win = 0.46*cos(pi*n) + 0.54;
  case 'Blackman'
    % Std Blackman: a = 0.16
    if nargin < 3, a=0.16; warning('Parameter "a" not supplied. "a" is set to 0.16.'); end %#ok
    n = 0.5+n/2; % rescale
    a0 = (1-a)/2; a1 = 1/2; a2 = a/2;
    win = a0-a1*cos(2*pi*n)+a2*cos(4*pi*n);
    %win = 1-win;
  case 'Tukey'
    if nargin < 3, error('parameter "a" not supplied. (win_name, n, a)'); end
    win = (n>=0).*((n<=(1-a/2)) +...
                   (n>1-a/2).*(0.5+0.5*cos(2*pi/a*(n-1+a/2)))...
                   ).*(n<1);
  otherwise
    error(sprintf('Uknown window function "%s".\nChoose: Tukey, Hanning, Hamming or Blackman.', win_name)) %#ok
end
end

