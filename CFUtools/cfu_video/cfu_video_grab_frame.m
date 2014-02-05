function cfu_video_grab_frame (settings, frame_idx_in)
%
% cfu_video_grab_frame(settings [,frame_idx])
% Grabs a video frame from figure and saves it as PNG.
%
% When frame_idx is given it forces the index of the current frame, when not given an
% internal counter is used.
%
% 2011-11-22, MFR, Version 1.0
% 2012-06-08, MFR, Version 1.1, Added frame_idx as input parameter.
% 2013-05-02, MFR, Version 1.2, Now uses saveSameSize to save video frames in the
%                               correct resolution.
% 2013-05-06, MFR, Version 1.3, frame_idx is now optional.
%

if nargin < 1, help(mfilename); return; end;

if nargin < 2, 
    % get the global variable
    global cfu_video_frame_idx
    % make sure it already existed (=not empty)
    if isempty(cfu_video_frame_idx)
        error(['The global variable ''cfu_video_frame_idx'' set by ' ...
               'cfu_video_init does not exist.'])
    else
        % Count up the frame idx
        cfu_video_frame_idx = cfu_video_frame_idx+1;
        frame_idx = cfu_video_frame_idx;
    end
else
    frame_idx = frame_idx_in;
end

file = sprintf('%s%s_frame%04i.png',settings.folder, settings.filename, frame_idx)

set(settings.hf,'Position',[1 1 settings.resolution]);
set(settings.hf,'Units','pixels');
drawnow
saveSameSize(settings.hf, 'format', 'png', 'file', file, 'dpi',settings.dpi);



