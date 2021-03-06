% Generate Spatial matched filter using slurm for parallel processing
% 

resolution = [5405,192];

adm.scriptname='SMF_slurm_generateline';
par.scanline = 49:144;%1:resolution(2);
% par.resolution = {resolution};
adm.slurm_opt = '--cpus-per-task=6 -xfcfu[10-12]';
cfu_cluster(adm, par);

% HOW TO
% 1. Setup run_cluster_function, script name and variables in par structure
% 2. Log in to fcfu1 - ssh fcfu1
% 3. Start Matlab and execute run_cluster_function
% Run_cluster_function will now start cfu_cluster and execute the script
% name suppled in the adm structure with parameters from par.
