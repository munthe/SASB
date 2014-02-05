% version = get_ffmpeg_version([quiet=0])
% Returns 'GID' or 'STD'
%

function version = get_ffmpeg_version(quiet)
if nargin < 1
    quiet = 0;
end
[ret_val return_str] = system('LD_LIBRARY_PATH=/usr/lib/; ffmpeg -version');
if ret_val ~= 0, 
    err_txt = sprintf(['You don''t seem to have ffmpeg installed. ' ...
                       'You can install ffmpeg by runing:\n' ...
                       ' sudo apt-get install ffmpeg']);
    version = -1;
    error(err_txt);
end

if ~isempty(findstr('git', return_str))
    version = 'GIT';
    if ~quiet
        fprintf('Using developer GIT FFmpeg.\n')
    end
elseif ~isempty(findstr('ubuntu0.12.04.1', return_str))
    version = '12.04';
    if ~quiet
        fprintf('Using Stock Ubuntu 12.04 FFmpeg.\n')
    end
elseif ~isempty(findstr('SVN', return_str))
    version = '10.04';
    if ~quiet
        fprintf('Using Stock Ubuntu 10.04 FFmpeg.\n')
    end
elseif ~isempty(findstr('lucid', return_str))
    version = '10.04';
    if ~quiet
        fprintf('Using Stock Ubuntu 10.04 FFmpeg.\n')
    end
    
else
    error('Could not determine FFmpeg version.')
end
end
