function st = cfu_tools_get_info()
% st = cfu_tools_get_info()
% Returns a struct containing info about the CFUtools repository.
%
% The struct contains the following members:
%      revision : Revision of the checked out repository
%          date : Date of the last (checked in) change.
%          path : path to CFUtools 
%           url : location of repository on the server
%   last_author : Last author to edit in the repository
%          uuid : universally unique identifier of the repository
%
% 2013-06-25, Mofi, Init version.
%

dir_ending = ['div' filesep 'cfu_tools_get_info.m'];
cfutools_dir = strtrim(which(mfilename));
cfutools_dir = cfutools_dir(1:end-length(dir_ending));

[val str] = system(sprintf('svn info %s', cfutools_dir));
if val ~=0, error('Could not get info on the CFUtools repository.'); end

str_cell = strread(str, '%s'); %string to cell array

st.revision = str_cell{12};
st.date     = [str_cell{29} '_' str_cell{30}]; 
st.path     = str_cell{2};
st.url      = str_cell{4};
st.last_author = str_cell{21};
st.uuid     = str_cell{10};

