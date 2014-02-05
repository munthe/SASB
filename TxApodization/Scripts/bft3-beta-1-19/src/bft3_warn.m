% function bft3_warn(varargin)
%
% print concise warning message
%
% $Id: bft3_warn.m,v 1.3 2011-03-10 18:43:17 jmh Exp $

% @file  bft3_warn.m
% @brief Warning function
% ======================================================================
function bft3_warn(varargin)
[name line] = bft3_caller_name(1);
fprintf(['Warn: %s %d: ' varargin{1}], name, line, varargin{2:end})
