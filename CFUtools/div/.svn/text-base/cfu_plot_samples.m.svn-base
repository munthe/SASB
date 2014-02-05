function cfu_plot_samples (root_dir, varargin)
% cfu_plot_samples (root_dir, varargin)
%
% Input: 
%  required: root_dir -path to the folder containing the subfolder 'seq_0001/'.
%  optional: key-value pair of: seq_idx, ems_idx, caxis.
%
% 2013-06-25, Mofi, Init Version.
%

opt.seq_idx = 1;
opt.ems_idx = 1;
opt.caxis   = [-50 50];
opt = cfu_parse_input_parameters(opt, varargin);


if nargin < 1, error('''path'' must be set.'); end
    


%% get list of all sub-folders
fprintf('Recursively searching the file system for measurements...\n')
folders = get_list_of_folders(root_dir);



%% get main folders that contain sequences
[folders num_seq] =  get_list_of_all_sequences(folders);
fprintf('Found %i measurements.\n\n', length(folders))

for idx = 1:length(folders)
    folder = folders{idx};
    num_emissions = get_num_ems(folder);
    if opt.seq_idx > num_seq(idx)
        warning(sprintf('sequence index set too high. Reducing seq_idx to %i', num_seq(idx)))
        seq_local = num_seq(idx);
    else
        seq_local = opt.seq_idx;
    end
        
    fprintf('Plotting emission: %3i/%i\n', opt.ems_idx, num_emissions);
    fprintf('Sequence:          %3i/%i\n', opt.seq_idx, seq_local);
    fprintf('Of measurement:\n%s\n\n', folder);
    
    % set the title equal the last 3 sub-folder names
    str_idx   = findstr(folder, filesep);
    if length(str_idx) > 2    
        start_idx = 1+str_idx(end-2);
    else
        start_idx = 1;
    end
    title_str = folder(start_idx:end);
    
    
    
    %% load and plot the data
    load(sprintf('%s/seq_%04i/elem_data_em%04i.mat', folder, opt.seq_idx, opt.ems_idx))

    smax   = max(samples, [], 2);
    smin   = min(samples, [], 2);

    figure;
    % plot the min and max values
    h1 = subplot(1,2,1);
    h = plot(smax, 1:size(samples,1));
    hold on
    h = plot(smin, 1:size(samples,1), 'r');
    xlim([-2^11-11 2^11+10])
    ylim([0 size(samples,1)])
    set(gca, 'XTick', [-2^11 2^11-1])
    set(gca, 'XTickLabel', ['-2^11 '; '2^11  '])
    ylabel('Sample number')
    set(get(h,'Parent'), 'YDIR','reverse')
    % make the min/max plot smaller
    pos1 = get(h1, 'position');
    pos1 = pos1 - [-0.01 -0.01 .2 0.01];
    set(h1, 'position',pos1)

    % plot the samples
    h2 = subplot(1,2,2);
    imagesc(samples, opt.caxis)
    % make a colorbar and move its label to save space
    cbh = colorbar;
    xlabel('Channel number')
    ylh = ylabel(cbh, 'Amplitude');
    ylpos = get(ylh,'position');
    set(ylh,'position',ylpos -[3 0 0])
    % remove the y-axis
    set(gca, 'YTick', [])
    %set the title
    title(title_str, 'interpreter','none')

    
    % move the samples-plot to the left and make it wider.
    pos2 = get(h2, 'position');
    pos2 = pos2 - [.25 -0.01 -.25 0.01];
    set(h2, 'position',pos2)

end


end  %END of FUNCTION




%%
%% Helper functions taken from CFU_MAKE_MEAN_SAMPLES
%%

function folders = get_list_of_folders(root_dir)
pathstr = genpath(root_dir);
seplocs = findstr(pathstr, pathsep);
first   = [1 seplocs(1:end-1)+1];
last    = seplocs(1:end)-1;
folders = arrayfun(@(a,b) pathstr(a:b), first, last, 'UniformOutput', false);
end



function [main_folder num_seq] =  get_list_of_all_sequences(folders)
%% Remove any existing "mean_seq" folder from the folders-list
main_folder = [];
num_seq     = [];

has_mean_seq = strfind(folders, 'mean_seq');
if ~isempty(has_mean_seq)
    folders_out = [];
    for folder_idx = 1:length(folders)
        if isempty(has_mean_seq{folder_idx})
            folders_out{end+1}      = folders{folder_idx};
        end
    end
    folders = folders_out;
end

% Find foders that contain sequence folders and count the number of sequences they contain
has_seq = strfind(folders, 'seq_');


main_folder_idx=1;
last_had_seq=0;
last_had_no_seq = 1;
for folder_idx = 1:length(folders)
    if isempty(has_seq{folder_idx}) % This folder contains no sequence folders
        folder_candidate = folders{folder_idx};
        if last_had_seq
            main_folder_idx = main_folder_idx+1;
        end
        num_seq(main_folder_idx) = 0;
        last_had_seq    = false;
        last_had_no_seq = true;
    else       % This folder does contains sequence folders
        if last_had_no_seq
            main_folder{main_folder_idx} = folder_candidate;
            num_seq(main_folder_idx) = 1;
        else
            num_seq(main_folder_idx) = num_seq(main_folder_idx)+1;
        end
        last_had_no_seq = false;
        last_had_seq    = true;
    end
end
end



function num_emissions = get_num_ems(root_dir)
% Get number of emissions
num_emissions = ls ([root_dir filesep 'seq_0001' filesep]); %get string of files
num_emissions = strread(num_emissions, '%s'); %string to cell array
num_emissions = size(num_emissions,1);        %size = number of files
end
