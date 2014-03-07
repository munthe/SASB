% Generate one line of the spatial matched filter, for use with SMF_slurm.m

addpath('../cluster');
addpath(genpath('../lib'));
addpath('./Scripts');
init_parameters;
savepath = '/data/cfudata3/mah/Spatial_matched_filter/';
tmpdir = 'tmp_SMF/';

line = par.line;
resolution = [5,3];

SMFline = Generate_SMF_line(line,resolution,useCaseParams,transducerType);

%% Saving SMF

if(isdir([savepath 'tmp_SMF']))
    rmdir([savepath 'tmp_SMF'],'s')
    pause(1)
end
mkdir([savepath tmpdir])

save([savepath tmpdir 'SMF_line_' num2str(line)],SMFline)
