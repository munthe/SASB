function cfu_cluster_run_matlab(varargin)
    disp(sprintf('MATLAB started at %s', datestr(now)))
    [s,r]=system('hostname');
    disp(sprintf('Running on host %s', r))
    clear s r
    for idx=4:2:nargin
        par.(varargin{idx}) = varargin{idx+1};
    end
    if varargin{3}==1
        feval(varargin{1}, par);
    else
        run(varargin{1})
    end
    