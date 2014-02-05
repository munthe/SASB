function cfu_video_encode(settings, clean_frames)
%
% Usage: cfu_video_encode(settings [,clean_frames=1])
%
% By Mofi
% Version 1.1, tor 02 maj 2013 10:37:03 CEST, Added test for which version of FFmpeg is being used.
% Version 1.2, fre 03 maj 2013 11:27:57 CEST, Frame rate can now be set with two decimals.
% Version 1.3, ons 08 maj 2013 11:51:59 CEST, Added loglevel switch.
%

if nargin < 2, clean_frames = 1; end
if nargin < 1, help(mfilename); return; end



%% Make sure ffmpeg is installed
FFmpeg_v = get_ffmpeg_version;
test_for_bash;

if strcmp(settings.encoding, 'H264')
    % make sure we can encode h.264
    retval = get_ffmpeg_h264_capability;
    if retval ~= 0, 
        err_txt = sprintf(['The installed FFmpeg does not have the h.264 ' ...
                           'codec. Try instead on fcfu1, or\n' ...
                           'see "FFmpeg" on the CFUWiki for how to compile FFmpeg with h264 support.']);
        error(err_txt);
    end
end

%% Make sure images are flatten before feeding them to FFmpeg
%fprintf('Flattening PNG images by using GhostScript (gs).\n')
cmd = sprintf(['for i in $(ls %s%s_frame*.png); do '...
               'convert $i -background black -flatten +matte $i; done'],... 
              settings.folder, ...
              settings.filename);
%[retval str] = system(cmd);
% $$$ if retval
% $$$     warning(['You do not have the command ''convert'' available on your system. ''convert'' ' ...
% $$$              'is provided by GhostScript (gs). The resulting movie may be larger than ' ...
% $$$              'necessary and (by some reports) have a worse quality.']);
% $$$ end




%% Set video codec
if strcmp(settings.encoding, 'MPEG4')
    encoder_options = '-vcodec mpeg4 ';
    
    file_ending = 'mp4'; %for increased portability use avi instead of mp4.
elseif strcmp(settings.encoding, 'MPEG1')
    encoder_options = '-vcodec mpeg1video ';
    file_ending = 'mpeg';
elseif strcmp(settings.encoding, 'H264')
    encoder_options = '-vcodec libx264 ';
    encoder_options = [encoder_options '-preset veryslow '];  %better compression 
    encoder_options = [encoder_options '-profile:v baseline '];  %better backward compatibility    
    file_ending = 'mp4';
end    
% Diasable audio
encoder_options = [encoder_options ' -an '];
% Set old type pixel format to be playable on older machines (and QuickTime)
encoder_options = [encoder_options ' -pix_fmt yuv420p '];
% Set scaling depending on which version of FFmpeg is used.
if strcmp(FFmpeg_v, 'GIT')
    scaling = '-q:v 4';
else % STD version
    scaling = '-qscale 4';
end
encoder_options = [encoder_options scaling];


% Tell FFmpeg where the libraries are located, if using standard
% Ubuntu version
library_link = '';
if isunix && ~strcmp(FFmpeg_v, 'GIT')
    library_link = 'LD_LIBRARY_PATH=/usr/lib/;';
end


%% set the log level
% Set loglevel flag depending on which version of FFmpeg is used.
if strcmp(FFmpeg_v, '10.04')
    loglevel = sprintf('-v %.0f', settings.loglevel);
else % git version 
    switch settings.loglevel
      case 0 
        loglevel = 'quiet ';    % Show nothing at all; be silent. 
      case 1 
        loglevel = 'fatal ';   % Only show fatal errors. 
      case 2 
        loglevel = 'error ';    % Show all errors, including ones which can be recovered from. 
      case 3 
        loglevel = 'warning '; % Show all warnings and errors.
      case 4
        loglevel = 'info ';    % Show informative messages during processing.
      otherwise error('loglevel must be an integer from 0 to 4.')
    end
    loglevel = ['-loglevel ' loglevel];
end


%% Metadata
comment = 'Made with CFU_video, by Morten F. Rasmussen';
metadata = sprintf(['-metadata title="%s" -metadata author="%s" -metadata publisher="%s" ' ...
                    '-metadata copyright="%s" -metadata comment="%s" -metadata year="%s" ' ...
                    '-metadata date="%s" -metadata artist="%s" '], ...
                   settings.title, settings.author, settings.publisher, settings.copyright, ...
                   comment, settings.year, settings.year, settings.author);


%% Encode
settings.folder = strrep(settings.folder, ' ', '\ '); % handle space in path
encode_cmd = sprintf(['%s ffmpeg -f image2 -r %.2f %s ' ...
                    '-i %s%s_frame%%04d.png ',...
                    '%s -y %s ', ...
                    '%s%s.%s\n'], ...
                     library_link, ...
                     settings.frame_rate, ...
                     loglevel, ...
                     settings.folder, ...
                     settings.filename, ...
                     encoder_options, ...
                     metadata, ...
                     settings.folder, ...
                     settings.filename, ...
                     file_ending);
%fprintf('Running: ''%s''\n', encode_cmd);
fprintf('Encoding with FFmpeg..\n')
[ret_val] = system(encode_cmd);

if ret_val
    err_str = sprintf(['\nError: an error occured when encoding the video.\nThe command ' ...
                       'called was:\n%s\n'], encode_cmd);
    if loglevel < 4
        err_str = [err_str 'If not enough information was printet, try with ''loglevel'' set to 4'];
    end 
    error(err_str);
end
if clean_frames
    fprintf('Deleting PNG files..\n');
    % delete images
    if isunix
        cmd2 = sprintf('rm %s%s_frame*.png ',settings.folder,settings.filename);
    else
        cmd2 = sprintf('delete %s%s_frame*.png ',settings.folder, settings.filename);
    end
    system(cmd2);
end

% Clear global variable
clear cfu_video_frame_idx

% close figure
try 
    close(settings.hf)
catch err
end
    

fprintf('Done!\n');
end













%% Helper functions




function test_for_bash
[ret_val str] = system('echo $SHELL');
if strcmp(str, sprintf('/bin/bash\n')) ~= 1
    warning ('You seem not to be using the Bash shell. ')
end
end

