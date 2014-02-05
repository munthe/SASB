%
% ret_val = get_ffmpeg_h264_capability
%
% Returns 0 when FFmpeg can encode h.264, not 0 otherwise.
%

function ret_val = get_ffmpeg_h264_capability
version = get_ffmpeg_version(1);
if strcmp(version, '10.04')
    codecs = '-formats ';
else
    codecs = '-codecs ';
end


% 'DEV' = Decoding, Encoding, Video -capability
[ret_val b] = system(sprintf(['LD_LIBRARY_PATH=/usr/lib/; ffmpeg %s 2>/dev/null | grep h264 | grep DEV'], ...
                             codecs));

end
