function st = cfu_video_init(filename, varargin)
%
% settings = cfu_video_init(filename, varargin)
%
% Input:
%   filename:   name of video file without file ending. (Must be set)
%   folder:     Full path to the directory where the video should be saved. 
%               If not set, it defaults to the folder 'video' under the current
%               directory.
%   resolution: Size of figures (and thereby th video) in
%               pixels. Defaults to [1024 768]. 
%   frame_rate: Frame rate of video. Defaults to 3.
%   visible:    Whether or not figures should be visible during frame grabbing. 
%               ('off'|'on') Defaults to 'on'. For this to work the following 
%               arguments must be passed 
%               to the plotting function:  "'parent', settings.hax".
%   encoding:   Which encoding to use. MPEG4 is standard. 
%               ('MPEG4'|'MPEG1'|'h264')
%   dpi:        Sets the dpi of the movie. Default is 0=screen resolution. Changing
%               this value overrules the meaning of the 'resolution' argument.
%   loglevel:   determines the amount of information logged (printet) by FFmpeg. 
%               Standard is 3=show errors+warnings. Accepted values are:
%               0=quiet, 1=fatal, 2=errors, 3=warnings, 4=info
%   title:      Metadata title of the video. If not set, it defaults to the filename.
%   author:     Metadata author of the video. If not set, it defaults to 'CFU'.
%   publisher:  Metadata publisher of the video. Usually does not need to be changed.
%   license:    Metadata license of the video. Usually does not need to be changed.
%   year:       Metadata year of the video. Usually does not need to be changed.
%
%
% Output:
%   settings: struct containing all settings.
%
% Example:
%    settings = video_init('video1', 'resolution', [1280 720], 'frame_rate', 20, ...
%                           'encoding', 'MPEG4', 'author', 'Morten Fischer Rasmussen',...
%                           'title', 'Cyst phantom');
%    Call 'edit cfu_video_example_mpeg4' to see how the cfu_video can be used.
%
%
% 2011-11-22, MFR, Version 1.0, Init version.
% 2012-03-13, MFR, Version 1.1, Force save folder to end with '/'.
% 2013-04-30, MFR, Version 1.2, Updated relative path test.
% 2013-05-01, MFR, Version 1.3, Added h.264 codec and made the script more Windows
%                               friendly (not tested).
% 2013-05-06, MFR, Version 1.4, Asks the user whether old movie frames in the output
%                               folder should be deleted.
% 2013-05-07, MFR, Version 1.5, Now accepts an optional settings 'dpi'. Setting 'dpi'
%                               to anything other than zero, results in the movie not
%                               having the resolution specified.
% 2013-05-08, MFR, Version 1.6, Now accepts the metadata flags: title, author, year,
%                               publisher and copyright. 
%

if nargin < 1, help(mfilename), return; end


st.folder      = [];
st.filename    = [];
st.resolution  = [];
st.frame_rate  = 10;
st.visible     = 'off';
st.encoding    = 'MPEG4';
st.dpi         = 0;
st.loglevel    = 3; %show only errors and warnings
st.title       = [];
st.author      = 'CFU';
st.publisher   = 'Center for Fast Ultrasound Imaging (CFU) at DTU, Denmark';
st.copyright   = ['This file may be shared with any and all, as long as credit to CFU is ' ...
                  'maintained.'];
st.year        = [];



%% ------------------------------------------------------------------
% Input Tests
% -------------------------------------------------------------------
% Verify and update input arguments
if nargin > 1
    st = cfu_parse_input_parameters(st,varargin);
end

% Resolution
if isempty(st.resolution)
    st.resolution = [1024 768];
else
    dim = size(st.resolution);
    if length (dim) ~= 2, 
        error ('resolution must be a 1x2 vector.');
    end
end
% format strings
st.encoding = upper(st.encoding);
st.visible  = lower(st.visible);
% test the resolution
if strcmp(st.encoding,'H264') && mod(st.resolution(2), 2)
    error(['The h.264 encoding ONLY accepts resolutions where the width is devisable by 2.'])
end

% test loglevel
if ~isnumeric(st.loglevel) || st.loglevel > 4 || st.loglevel < 0
    error('loglevel must be an integer between 0 and 4.')
end

% Test the frame rate
if strcmp(st.encoding,'MPEG1') && st.frame_rate < 20
    error(['The minimum frame rate of MPEG1 is 20 fps.'])
end

% Set Folder
if isempty(st.folder)
    st.folder = [pwd filesep 'video' filesep];
elseif ~strcmp(st.folder(1),filesep) && ...
        ~strcmp(st.folder(1:2), '~/') && ...
        ~strcmp(st.folder(1:2), 'C:')  
    % Then it must be a relative path
    % Make it to an absolute path
    st.folder = [st.folder]; 
end
% handle empty space in path
%st.folder = strrep(st.folder, ' ', '\ ');

% set title
if isempty(st.title)
    st.title = filename;
end


% force ending with '/'
if st.folder(end) ~= filesep
    st.folder(end+1) = filesep;
end

%make folder, if it does not exist
if ~exist(st.folder, 'dir')
    retval = mkdir(st.folder);
    if retval ~= 1
        error(sprintf(' Could not create folder: %s\n', st.folder))
    end
end



% encoding
if ~strcmp(st.encoding, 'MPEG1') && ...
        ~strcmp(st.encoding, 'MPEG4') && ...
        ~strcmp(st.encoding, 'H264')
    error(' Error: encoding must be set to either ''MPEG1'', ''MPEG4''', ...
          ' or ''H264''.')
end


%framerate test
if strcmp(st.encoding, 'MPEG1') &&  st.frame_rate < 20
    error('MPEG1 does not support frame rates lower than 20 fps.')
end


% Meta data -- title
if isempty(st.title)
    st.title = settings.filename;
end
% Meta data -- date
if isempty(st.year)
    date_str = date;
    year = date_str(end-3:end);
    st.year = year;
end



% filename
st.filename = filename;



%% ------------------------------------------------------------------
% Create Figure and set its size
% -------------------------------------------------------------------
st.hf = figure('visible', st.visible);
st.hax = axes;
set(st.hf,'Position',[1 1 st.resolution]);
drawnow;




%% ------------------------------------------------------------------
% Clean old movie frames
% -------------------------------------------------------------------
% Make sure there doesn't exist old movie frames. If there does, ask the user whether 
% they should be deleted.
if old_frames_exists(st)
    prompt = ['There already exists old movie frames in the output directory.\n Do you want ' ...
              'to delete them? (Y/n)\n'];
    str = input(prompt,'s');
    str = lower(str);
    % Set 'y' as default
    if isempty(str)
        str = 'y';
    end
    if strcmp(str, 'y')
        clear_old_frames(st);
        fprintf('Deleted old movie frames.\n')
    elseif strcmp(str, 'n')
        fprintf('Leaving old movie frames untouched.\n')
    else
        error(sprintf('You must type ''y'' or ''n''.\n'));
    end
end





%% ------------------------------------------------------------------
% frame_idx
% -------------------------------------------------------------------
% make frame index a global variable
clear cfu_video_frame_idx
global cfu_video_frame_idx
cfu_video_frame_idx = 0;



end  %END of cfu_video_init











%% ------------------------------------------------------------------
% Helper Functions
% -------------------------------------------------------------------

function exists = old_frames_exists(st)
% DETECT_SAVED_FRAMES - Tests for already saved frames
%   Returns 1 when frames are detected, 0 otherwise.

if isunix
    ls = 'ls';
else
    ls = 'dir';
end
cmd = sprintf('%s %s%s*.png', ls, st.folder, st.filename);
[retval str] = system(cmd);

exists = ~retval;
end





function clear_old_frames(st)
% CLEAR_OLD_FRAMES - Deletes old frames.
%   
if isunix
    rm = 'rm -rf';
else
    rm = 'rem';
end
cmd = sprintf('%s %s%s*.png', rm, st.folder, st.filename);
[retval str] = system(cmd);
end
