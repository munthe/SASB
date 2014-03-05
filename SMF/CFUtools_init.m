%% Add projects to the path
% Path to CFUtools
cfupath = '../CFUtools';

% mfr_emission
addpath([cfupath filesep])

% Onda Library
%addpath([cfupath filesep 'onda_lib' filesep 'src'])
addpath([cfupath filesep 'onda_lib2'])

% cfu_video -- create MPEG1, MPEG4 and MPEG4+h.264 Videos from within MATLAB
addpath([cfupath filesep 'cfu_video'])

% matlabfrag -- Print LaTeX and EPS graphics
addpath([cfupath filesep 'matlabfrag'])

% mlf2pdf -- matlabfrag-to-PDF -- Print LaTeX and PDF graphics
addpath([cfupath filesep 'mlf2pdf'])

% real2rgb -- linear colorscales in bot RGB and BW
addpath([cfupath filesep 'real2rgb'])

% export_fig -- This function saves a figure or single axes to one or more vector and/or bitmap file formats, and/or outputs a rasterized version to the workspace
addpath([cfupath filesep 'export_fig'])

% Load Field II
os = computer;
switch os
  case {'PCWIN' , 'PCWIN64'}
    addpath([cfupath filesep 'Field_win']);
  case {'GLNX86', 'GLNXA64'}
    addpath([cfupath filesep 'Field_linux']);
  case 'MACI64'
    addpath([cfupath filesep 'Field_macosx']);
  otherwise
    error('CFUtools error: Your operating system was not recognized.')
end
% Additional  (small) scripts
addpath([cfupath filesep 'div'])

% Rubberband -- used to select regions in an image with the mouse. Used by
% cfu_get_psf_metrics and cfu_get_cyst_metrics.
addpath([cfupath filesep 'rubberband'])


% freeColors
addpath([cfupath filesep 'freezeColors'])
