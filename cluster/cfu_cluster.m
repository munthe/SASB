% Run a number of jobs automatically. This script is intended for use
% on a cluster running SLURM, but it can also be used in a single
% MATLAB instance if necessary.
%
% Note that all numbers are converted to strings and back, as they
% need to go through the shell to get into SLURM. Therefore, avoid
% decimal numbers as far as possible. A smarter implementation may
% be possible, but has not yet been investigated.
%
% The function has two input arguments:
%
% - An administration struct with the following fields:
%
%   scriptname - (mandatory) string, name of the MATLAB
%   script/function (user script) to run on the cluster.
%
%   is_function - (optional) logical, true if the user script is a
%   function, otherwise false. If the user script is a function,
%   see function_args. Otherwise, see script_args. Default false
%
%   function_args - (optional) string, specify how the user script
%   accepts arguments, when the user script is a function.
%     'struct' - The user script takes a single argument, which is
%     a struct containing the parameters
%     'named' - The user script takes named parameters, e.g.,
%     my_func('c', 1480, 'fs', 100e6) - NOT IMPLEMENTED
%     'list' - The user script takes its arguments in a specified
%     order, e.g., my_func(1480, 100e6) - NOT IMPLEMENTED
%   Default 'struct'
%
%   script_args - (optional) string, specify how the user script
%   accepts arguments, when the user script is not a function.
%     'struct' - The workspace will contain a structure with the
%     name 'par' containing all the parameters
%     'nostruct' - The workspace will contain a variable for each
%     parameters. The name of the variables are the names of the
%     fields in the parameter struct (see below). - NOT IMPLEMENTED
%   Default 'struct'
%
%   slurm - (optional) logical, run using SLURM (in parallel) or in
%   the calling MATLAB session (sequential), default true (use
%   SLURM). If true, the MATLAB session must run on one of the cfu
%   or fcfu machines.
%
%   slurm_opt - (optional) string, options passed to sbatch,
%   default ''
%
%   save_params - (optional) logical, save a mat file with an array of
%   structs containing all combinations of parameters the user script
%   is called with, default false
%
%   save_params_filename - (optional) string, name of mat file to save
%   the run parameters to, if save_params is true, default
%   cfu_cluster.mat
%
%   no_run - (optional) logical, do not actually run the user
%   script. This can be used for debug, e.g., in combination with
%   save_params, default false
%
%   is_slurm - (optional) logical, put a field is_slurm in the
%   parameter struct for the user script. The value of this field is
%   true if the user script runs under SLURM, otherwise false. It can
%   be used to terminate MATLAB in SLURM when execution has finished.
%
% - A struct or array of structs with the parameters to run the
%   specified script with. If a single struct is given, each field can
%   be a scalar, an array, or a cell array. The cartesian product is
%   then taken over all fields to create the sets of parameters to run
%   the user script with. If an array of structs is given, each struct
%   is assumed to contain the parameters for a single instance of the
%   user script.
%
% The following example runs 'my_field_sim' with all combinations of
% the specified focal depths, transmit apodizations, and number of
% active elements.
%
%  adm.scriptname='my_field_sim';
%  par.r_focus=[20 40]/1000; % Depth of focus point
%  par.tx_apod={'rectwin', 'hamming'}; % TX apodizations
%  par.n_active=[32 64 128]; % Number of active elements
%  cfu_cluster(adm, par);
%
% The next example shows how to let the user script derive all the
% parameters. The only field in the par struct is an index, which
% the user script can then use to calculate the parameter set to
% run. This can be used to avoid issues when converting to and from
% strings, as is currently required.
%
%  adm.scriptname='my_field_sim';
%  par.index=1:128; % my_field_sim derives all parameters from this
%  cfu_cluster(adm, par);
%
% Version info:
%   2013/06/09, mbs, Initial version

function cfu_cluster(adm_user, par_user)
    
    if nargin~=2
        error('cfu_cluster takes exactly two arguments')
    end
    
    if ~isstruct(adm_user) || ~isstruct(par_user)
        error('Both arguments to cfu_cluster must be structs')
    end
    
    if ~isfield(adm_user, 'scriptname')
        error('Mandatory field ''scriptname'' missing from adm_user struct')
    end
    
    % Default values
    [status,hostname]=system('hostname');
    adm.scriptname='';
    adm.is_function=false;
    adm.function_args='struct';
    adm.script_args='struct';
    adm.slurm=true;
    if hostname(1:3)=='cfu'
        adm.slurm_opt='-x cfu[25] ';
    else
        adm.slurm_opt='';
    end
    adm.save_params=false;
    adm.save_params_filename='cfu_cluster.mat';
    adm.no_run=false;
    adm.is_slurm=true;
    
    % Names of fields in adm_user struct
    adm_user_names=fieldnames(adm_user);
    
    % Update adm struct with user values
    for idx=1:length(adm_user_names)
        if ~isfield(adm, adm_user_names{idx})
            error(['Unsupported field in adm struct: ' adm_user_names{idx}])
        end
        if strcmpi(adm_user_names{idx}, 'slurm_opt')
            adm.(adm_user_names{idx}) = [adm.(adm_user_names{idx}) adm_user.(adm_user_names{idx})];
        else
            adm.(adm_user_names{idx}) = adm_user.(adm_user_names{idx});
        end
    end
    
    % A few things still have to be implemented
    if strcmpi(adm.function_args, 'named')
        error('Named function arguments are not yet implemented, option adm.function_args=''named''')
    elseif strcmpi(adm.function_args, 'list')
        error('List of function arguments are not yet implemented, option adm.function_args=''list''')
    elseif strcmpi(adm.script_args, 'nostruct')
        error('Individual variables in workspace are not yet implemented, option adm.script_args=''nostruct''')
    end
    
    % Create parameter struct array
    if length(par_user)==1
        % Get fields
        par_fields=fieldnames(par_user);
        
        % Get number of values for each field
        n_entries=zeros(1, length(par_fields));
        for idx=1:length(par_fields)
            if isstr(par_user.(par_fields{idx}))
                n_entries(idx)=1;
            else
                n_entries(idx)=length(par_user.(par_fields{idx}));
            end
        end
        
        % Get the total number of parameter combinations
        n_params=prod(n_entries);
        
        % Array of indexes used to iterate over
        idxs=ones(1, length(par_fields));

        % Create all parameter sets
        for p_idx=1:n_params
            for idx=1:length(par_fields)
                if iscell(par_user.(par_fields{idx}))
                    par(p_idx).(par_fields{idx})=par_user.(par_fields{idx}){idxs(idx)};
                else
                    if isstr(par_user.(par_fields{idx}))
                        par(p_idx).(par_fields{idx})=par_user.(par_fields{idx});
                    else
                        par(p_idx).(par_fields{idx})=par_user.(par_fields{idx})(idxs(idx));
                    end
                end
                if adm.is_slurm
                    par(p_idx).is_slurm=adm.slurm;
                end
            end
            
            % Increment idxs, principle of a carry-ripple adder
            carry=1;
            for idx=length(par_fields):-1:1
                if carry==1
                    if idxs(idx)<n_entries(idx)
                        idxs(idx)=idxs(idx)+1;
                        carry=0;
                    else
                        idxs(idx)=1;
                    end
                end
            end
        end
    else
        par=par_user;
    end
    par_all=par;
    clear par;
    
    % Save parameters
    if adm.save_params
        save(adm.save_params_filename, 'par_all');
    end
    
    % Should the user script be run?
    if adm.no_run
        return
    end
    
    % Should SLURM be used?
    if adm.slurm
        par_fields=fieldnames(par_all(1));
        for idx=1:length(par_all)
            tic
            slurm_str=['sbatch -N1 -n1 ' adm.slurm_opt ' cfu_cluster_srun_matlab.sh ' adm.scriptname ' is_function '];
            if adm.is_function
                slurm_str=[slurm_str '1 '];
            else
                slurm_str=[slurm_str '0 '];
            end
            for pidx=1:length(par_fields)
                slurm_str=[slurm_str par_fields{pidx} ' '];
                if isstr(par_all(idx).(par_fields{pidx}))
                    if strcmpi(par_all(idx).(par_fields{pidx}), '')
                        error('Empty strings are not supported for SLURM processing')
                    end
                    slurm_str=[slurm_str par_all(idx).(par_fields{pidx}) ' '];
                elseif isempty(par_all(idx).(par_fields{pidx}))
                    error('Empty variables are not supported for SLURM processing')
                else
                    slurm_str=[slurm_str num2str(par_all(idx).(par_fields{pidx})) ' '];
                end
            end
            str=sprintf('Submitting parameter set %d of %d', idx, length(par_all));
            disp(str);
            toc
            system(slurm_str);
            toc
        end
        disp('All SLURM jobs have now been submitted.');
        disp('You can check the status using the command squeue in a terminal');
        disp('on the appropriate cluster.');
    else
        for idx=1:length(par_all)
            str=sprintf('Running parameter set %d of %d', idx, length(par_all));
            disp(str);
            par=par_all(idx);
            if adm.is_function
                feval(adm.scriptname, par);
            else
                cfu_cluster_run_script(adm.scriptname, par);
            end
        end
    end
    