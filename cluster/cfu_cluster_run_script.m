%% Function to wrap running a script to have a clean workspace

function cfu_cluster_run_script(scriptname, par)
    disp(sprintf('MATLAB started at %s', datestr(now)))
    [s,r]=system('hostname');
    disp(sprintf('Running on host %s', r))
    clear s r
    run(scriptname);
    