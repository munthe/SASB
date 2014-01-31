% function [name, line] = bft3_caller_name(level)
%
% return name (and line) of calling routine or file (if level=1, the default)
% or name further up or down the stack by changing caller
%
% $Id: bft3_caller_name.m,v 1.5 2011-07-25 15:07:57 jmh Exp $

%> @file  bft3_caller_name.m
%> @brief Return name and line of calling routine
% ======================================================================
%> function [name, line] = bft3_caller_name(level)
%>
%> return name (and line) of calling routine or file (if level=1, the default)
%> or name further up or down the stack by changing caller
%>
function [name, line] = bft3_caller_name(level)

if nargin < 1
  level = 1;
end

st = dbstack;
if length(st) < 2+level
  name = '';
  line = 0;
else
  st = st(2+level);
  name = st.name;
  line = st.line;
end
