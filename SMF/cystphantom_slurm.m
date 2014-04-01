
adm.scriptname='cystphantom';
par.scanline = 1;
adm.slurm_opt = '--cpus-per-task=4 -xfcfu[10-12]';
cfu_cluster(adm, par);

% HOW TO
% 1. Setup run_cluster_function, script name and variables in par structure
% 2. Log in to fcfu1 - ssh fcfu1
% 3. Start Matlab and execute run_cluster_function
% Run_cluster_function will now start cfu_cluster and execute the script
% name suppled in the adm structure with parameters from par.
