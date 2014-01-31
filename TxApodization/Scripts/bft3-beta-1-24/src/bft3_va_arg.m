% @file  bft3_va_arg.m
% @brief function for updating a struct using (name, value) pairs
% ======================================================================
%
% ======================================================================
%> @brief function for updating a struct using (name, value) pairs
%>
%> @param opt
%> @param varargs
%>
%> @return opt
% ======================================================================
function opt = bft3_va_arg(opt,varargs)
%function opt = bft3_va_arg(varargin)
% splits argument varargin and returns a struct with (name,value) pairs,
% e.g. opt = bft3_va_arg(opt,{'a',2,'b',3}) returns
% opt =
%     a: 2
%     b: 3
%
% $Id: bft3_va_arg.m,v 1.4 2011-07-25 02:57:11 jmh Exp $
%
  if nargin < 1, help(mfilename), error(mfilename), end
  if nargin == 1 && strcmp(opt, 'test'), bft3_va_arg_test, clear, return, end

  % Allow introduction of new members in output (default: no)
  allow_new_members = 0;

  % Member names
  names = fieldnames(opt);

  % Un-even number of elements in varargs, last element is assumed to
  % be a bolean indicating whether we allow new members
  if (mod(length(varargs),2)==1)
    allow_new_members = varargs{end};
  end
  
  for k=1:1:floor(length(varargs)/2)
    eval(['arg','=varargs{',int2str(2*k-1),'};']);
    if ismember(arg,names)
      eval(['opt.',arg,'=varargs{',int2str(2*k),'};']);
      eval(['b = opt.',arg,';']);
      if isstruct(b)
        eval(['opt.',arg,' = bft3_update_struct(opt.',arg,[',' ...
                            'varargs{'],int2str(2*k),'});']);
      end
    else
      if allow_new_members
        disp(['New struct member ' arg ' introduced']);
        eval(['opt.',arg,'=varargs{',int2str(2*k),'};']);      
      else
        error(['Unknown struct member ' arg]);
      end
    end
  end
end

function bft3_va_arg_test
  args = {'a', 1, 'b', 2, 'c', 3};
  opt.a = 0;
  opt.b = 0;
  opt.c = 0;
  opt = bft3_va_arg(opt, args);
  
  if ((opt.a ~= cell2mat(args(2))) || (opt.b ~= cell2mat(args(4))) || (opt.c ~= cell2mat(args(6))))
    error 'Should not get here'
  end
end
