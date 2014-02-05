function opt = cfu_parse_input_parameters(opt,varargs)
%
% param_struct = CFU_PARSE_INPUT_PARAMETERS (param_std, param1_name, param1_value, ... )
%  or
% param_struct = CFU_PARSE_INPUT_PARAMETERS (param_std, varargs)
%
% Parses the input parameters given in parameter-value pairs. Any number of
% parameter-value pairs can be given.
%  
% Input:
%   param_std:    Struct containing the standard parameter values. Members can also be empty. 
%                 Only parameter/members defined in this struct can be updated.
%
%   param_name:   String containing parameter name to be updated. This name must match
%                 an existing 
%   param_value:  New value for the parameter name.
%
% Output:
%   para_struct:  Struct containing the updated values.
%
% Author: Morten Fischer Rasmussen
% 2013-05-08, MFR, Init version.
% 2013-08-29, MFR, Optimized for speed.
%

if mod(length(varargs),2) ~= 0
    error('parse_input_parameters:InvalidOptionalParameters', ...
          'Optional parameters must be specified as parameter-value pairs');
end

% Existing Member Names
names_std = fieldnames(opt);

% Input: member names
names_new  = {varargs{1:2:end}};
% Input: new values
values_new = {varargs{2:2:end}};

%make sure member exists
[member_exists member_idx] = ismember(names_new, names_std);

% Update values 
for idx = 1:length(names_new)
    if member_exists(idx)
        opt.(names_std{member_idx(idx)}) = values_new{idx};
    else
        err_str = sprintf('Parameter name ''%s'' does not exist.', names_new{idx});
        error(err_str)
    end
end
end
