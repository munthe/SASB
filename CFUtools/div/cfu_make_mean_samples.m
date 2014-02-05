function [folders_out num_seq_out]= cfu_make_mean_samples (root_dir)
% make_mean_samples (root_dir)
%
% Input: root_dir - path to root folder which contains the sequence folders 
%                  (e.g., 'root_dir/seq_0001/' exists.)
%
% The average of all sequences is saved in 'root_dir/mean_seq/'.
%
% 2013-06-25, Mofi, Init Version.
% 2013-09-08, Mofi, Version 1.1, Now recursively searches the sub-folders for
%                                measurements. 
% 2013-10-24, Mofi. Version 1.2, Now prints the current working directory.
%


fprintf('Recursively searching the file system for measurements...\n')

%% get list of all sub-folders
folders = get_list_of_folders(root_dir);


%% get main folders that contain sequences
[folders num_seq] =  get_list_of_all_sequences(folders);


%% remove folders that already contains a the mean sequencels
unfinished_avg = [];
remove_folder  = [];
for idx = 1:length(folders)
    if contains_mean_seq(folders{idx})
        num_ems_seq_mean = size(strread(ls([folders{idx} filesep 'mean_seq/']), '%s'),1);
        num_ems_seq_0001 = size(strread(ls([folders{idx} filesep 'seq_0001/']), '%s'),1);
        
        if num_ems_seq_mean == num_ems_seq_0001        
            remove_folder(end+1) = idx;
        else
            unfinished_avg(end+1) = idx;
        end
    end
end
if ~isempty(unfinished_avg)
    fprintf('Adding the following unfinished ');
    if length(unfinished_avg) > 1
        fprintf('averages to the working pool:\n');
        for idx =1:length(unfinished_avg)
            fprintf('%s\n', folders{idx});
        end
    else
        fprintf('average to the working pool:\n');
        fprintf('%s\n', folders{1});
    end
end
fprintf('\n')
%The actual removal
num_seq(remove_folder) = [];
folders(remove_folder) = [];


%% remove folders that only contain a single sequence
remove_folder = [];
for idx = 1:length(folders)
    if num_seq(idx) == 1;
        remove_folder(end+1) = idx;
    end
end
num_seq(remove_folder) = [];
folders(remove_folder) = [];


%% Print INFO
if length(folders) == 0
    fprintf('No measurements found.\n\n');
    return
else
    fprintf('Found %i measurements to average:\n',length(folders))
    for idx = 1:length(folders)
        fprintf('%3i seq: %s\n', num_seq(idx), folders{idx});
    end
    fprintf('\n');
end



%% Calculate and save the average
for avg_idx = 1:length(folders)
    num_ems = get_num_ems(folders{avg_idx});
    make_average(folders{avg_idx}, num_seq(avg_idx), num_ems);
end



% set output vars if needed;
if nargout ==2,
    folders_out = folders;
    num_seq_out = num_seq;
end
end









%%%%%--------------------------------------
% Helper/ Worker scripts




function make_average(root_dir, num_sequences, num_emissions)
% create the output dir
out_dir = [root_dir filesep 'mean_seq'];
cfu_mkdir(out_dir)

fprintf ('\nWorking on: %s\n', root_dir);
fprintf ('Number of emissions: %i\n', num_emissions);
fprintf ('Number of sequences: %i\n', num_sequences);


disp 'Emission: '
for ems_idx=1:num_emissions
    fprintf('%i ', ems_idx)
    for seq_idx=1:num_sequences
        load (sprintf('%s%sseq_%04i%selem_data_em%04i', ...
                      root_dir, filesep, ...
                      seq_idx,  filesep, ...
                      ems_idx))
        if ems_idx == 1 && seq_idx == 1
            compl_seq = zeros([num_sequences size(samples)]);
        end
        compl_seq(seq_idx, :,:) = double(samples);
    end
    
    % calc mean and save
    mean_samples = squeeze(mean(compl_seq));
    save(sprintf('%s%selem_data_em%04i.mat', out_dir, filesep, ems_idx), 'mean_samples')
end
fprintf('\n')
end





function num_emissions = get_num_ems(root_dir)
% Get number of emissions
num_emissions = ls ([root_dir filesep 'seq_0001' filesep]); %get string of files
num_emissions = strread(num_emissions, '%s'); %string to cell array
num_emissions = size(num_emissions,1);        %size = number of files
end


function num_sequences = get_num_seq(root_dir)
% Get number of sequences
num_sequences = ls (root_dir);                %get string of files and folders
num_sequences = strread(num_sequences, '%s'); %string to cell array
num_sequences = strfind(num_sequences, 'seq_');%gives 1 when found, empty when not. 
num_sequences = sum([num_sequences{:}]);      %cell array to vector,count number of ones
end


function flag = contains_mean_seq(root_dir)
% Get number of sequences
str = ls (root_dir);                %get string of files and folders
str = strread(str, '%s'); %string to cell array
str = strfind(str, 'mean_seq');%gives 1 when found, empty when not. 
flag = sum([str{:}]);      %cell array to vector,count number of ones
end





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