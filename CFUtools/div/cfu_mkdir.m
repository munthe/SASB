%
% function that creates a directory, and throws an error if not successfull.
% cfu_mkdir (dir)
%
% By Morten F. Rasmussen
% Version 1.0 Init version
%

function cfu_mkdir (dir)

if ~ischar(dir)
    error(' Error: input must be a string.');
end

if (exist(dir, 'dir') ~= 7)
    retval = mkdir (dir); 
    if (retval ~= 1) 
        error(sprintf('Could not create folder \"%s\"', dir));
    end
end
