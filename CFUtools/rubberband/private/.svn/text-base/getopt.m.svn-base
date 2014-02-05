% getopt -- parse command options.
%
%      [opt, arg] = getopt(spec, ...)
%
% The 'getopt' command parses options from an argument list.
%
% First argument SPEC is a structure array where each element describes
% an option.  An option structure has these fields:
%
% name
%      The name of the option (a string) preceded by a hyphen.
% type
%      The type of the option (a string).  This field determines whether
%      the option requires an argument or not.  The option type is either
%      one of the built-in option types (see below), a Matlab class name,
%      or a user-defined type.
% index
%      The structure array index of the option.  This field is only
%      mandatory for alias options.
% predicate
%      The predicate function for a user-defined option type.  Value is
%      either a function name, a function handle, or a cell array whose
%      first element is a function and the remaining elements are any
%      number of additional arguments for the predicate function.  When
%      calling the predicate function, the option argument to be checked
%      is passed as the first argument to the predicate function.  First
%      return value is the predicate.  Optional second return value is
%      the actual option value.
% choices
% values
%      The choices field contains a set of valid option arguments.  If the
%      choices field is empty, any argument of the correct type is accepted.
%      If the choices *and* values field is not empty, the option argument
%      is replaced by the corresponding choice value.
% value
%      The value of the option.   For options requiring an argument, the
%      option value is set to the option argument.
% assign
%      The variable name in the caller's workspace where the option value
%      will be assigned to.  If the value of this field is a valid Matlab
%      variable name, this variable name will be used as is.  Otherwise,
%      the variable name will be derived from the option name by replacing
%      all non-alphanumeric characters with underscores (the initial hyphen
%      is ignored).  For the later case, the variable is assigned iff the
%      value of the 'assign' field is true.
%
% The 'getopt' command returns the modified option structure array OPT,
% and the remaining arguments ARG.  Options and arguments are permuted while
% being parsed.  The special argument '--' terminates parsing in all cases.
%
% An option argument can be provided in two ways.  Either as a separate
% argument immediately following the option it belongs to, or together with
% the option separated from it by an equal sign.
%
% If the argument SPEC is empty or a number, 'getopt' does not accept any
% options and the return value OPT is an option structure array with that
% many elements.
%
% Built-in option types:
%
% alias
%      The option is an alternative name for another option.  The index
%      field is the structure array index of the original option.  If the
%      value field is not empty, an alias option does not take an argument
%      but sets the value of the original option to that value.
% flag
%      The option is a simple "has seen" option.  The option does not take
%      an argument and the option value defaults to zero.
% counter
%      The option is a counter and it does not take an argument.  Each time
%      the option occurs the value is incremented by one.  The option value
%      defaults to zero.
% toggle
%      The option is an on/off option.  Each time the option occurs, the
%      value is logically negated.  The option value defaults to zero.
%
% The following built-in option types can be overwritten by user-defined
% option types:
%
% color
%      The option requires a Matlab color specification as its argument.
%      Option argument can be a short or long Matlab color name, a RGB
%      intensity tripple, or a gray-scale intensity value.
% cursor
%      The option requires a Matlab cursor name as its argument.
% line
%      The option requires a Matlab line style as its argument.
%      An argument of 'solid', 'dashed', 'dotted', and 'dashdot' is
%      accepted as an alternative name for the corresponding Matlab
%      line style.
% marker
%      The option requires a Matlab marker symbol as its argument.
%      An argument of 'plus', 'circle', 'asterisk', 'star', 'point',
%      'cross', 'x-mark', 'square', 'diamond', 'up', 'down', 'left',
%      'right', 'pentagram', and 'hexagram' is accepted as an alternative
%      name for the corresponding Matlab marker symbol.
% symbol
%      The option argument must be a Matlab symbol name.
% number
%      The option argument must be a number.
% positive_number
%      The option argument must be a positive number, i.e. greater than
%      zero.
% negative_number
%      The option argument must be a negative number, i.e. less than zero.
% non_positive_number
%      The option argument must be a non-positive number, i.e. less than
%      or equal to zero.
% non_negative_number
%      The option argument must be a non-negative number, i.e. greater
%      than or equal to zero.
% non_zero_number
%      The option argument must be a non-zero number.
% integer
%      The option argument must be an integral number.
% positive_integer
%      The option argument must be a positive integral number.
% negative_integer
%      The option argument must be a negative integral number.
% non_positive_integer
%      The option argument must be a non-positive integral number.
% non_negative_integer
%      The option argument must be a non-negative integral number.
% non_zero_integer
%      The option argument must be a non-zero integral number.
% any
%      The option argument can be of any type.
%
% Example:
%
%      function foo(varargin)
%
%        % Initialize options.
%        opt = getopt(5);
%
%        for ind = 1:numel(opt)
%          switch ind
%           case 1
%            opt(ind).name = '-flag';
%            opt(ind).type = 'flag';
%           case 2
%            % Variable names equal to Matlab commands have to be
%            % initialized (work around Matlab 'assignin' bug).
%            mode = 4;
%            opt(ind).name = '-mode';
%            opt(ind).type = 'char';
%            opt(ind).choices = {'box', 'line'};
%            opt(ind).values = [4, 2];
%            opt(ind).value = mode;
%            opt(ind).assign = 1;
%           case 3
%            opt(ind).name = '-color';
%            opt(ind).type = 'color';
%           case 4
%            opt(ind).name = '-cursor';
%            opt(ind).type = 'cursor';
%           case 5
%            opt(ind).name = '-fullcross';
%            opt(ind).type = 'alias';
%            opt(ind).index = 4;
%            opt(ind).value = 'fullcrosshair';
%           otherwise
%            error('Should not happen');
%          end
%        end
%
%        % Parse options from argument list.
%        [opt, arg] = getopt(opt, varargin{:});

%% getopt.m --- the preceding comment is the documentation string.

% Copyright (C) 2004--2009 Ralph Schleicher

% Author: Ralph Schleicher <rs@ralph-schleicher.de>
% Time-stamp: <2009-06-20 16:33:47 CEST>

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License as
% published by the Free Software Foundation; either version 2,
% or (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program; see the file COPYING.  If not, write to
% the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA 02111-1307, USA.

% As a special exception, Ralph Schleicher gives permission to link
% the code of this program with MATLAB from The Mathworks, Inc. (or
% with modified versions of MATLAB that use the same license as
% MATLAB), and distribute linked combinations including the two.
% You must obey the GNU General Public License in all respects for
% all of the code used other than with MATLAB.  If you modify this
% file, you may extend this exception to your version of the file,
% but you are not obligated to do so.  If you do not wish to do so,
% delete this exception statement from your version.

%% Commentary:

% This code is written with portability in mind.  Therefore, all
% M-Lint warnings are errors of M-Lint in the first place, not mine.
% And yes, I refuse cluttering my code with '%#ok' labels to make
% M-Lint happy.

% Members of the Church of Emacs may appreciate the following code
% snippet.
%
% (defun matlab-renumber-switch (start end)
%   "Renumber the labels of a Matlab switch statement."
%   (interactive "r")
%   (save-excursion
%     (goto-char start)
%     (let ((i (prefix-numeric-value current-prefix-arg)))
%       (while (re-search-forward "^\\(\\s-*case\\)\\s-*[0-9]*\\s-*$" end t)
%         (replace-match (format "\\1 %d" i) t)
%         (setq i (1+ i))))))

%% Code:

% Program entry point.
function [opt, rest] = getopt(spec, varargin)

  % Check number of arguments.
  error(nargchk(1, inf, nargin));

  % Non-empty means to create an option structure array
  % with that many elements.
  elem = [];

  if isstruct(spec)
    0;
  elseif isempty(spec)
    elem = 0;
  elseif isnumeric(spec)
    if ~ getopt_non_negative_integer_p(spec)
      error('Option structure size has to be a non-negative integral number');
    end
    elem = spec;
  else
    error('Invalid argument');
  end

  if ~ isempty(elem)
    if elem == 0
      c = cell(0);
    else
      c = cell(elem, 1);
    end
    % Create option structure.
    opt = struct('name', c, ...
                 'type', [], ...
                 'index', [], ...
                 'predicate', [], ...
                 'choices', [], ...
                 'values', [], ...
                 'value', [], ...
                 'assign', []);
    % Return remaining arguments as is.
    rest = varargin;
    return;
  end

  % Check option structure array.
  opt = spec;

  % Add missing fields.
  fields = fieldnames(opt);
  switch 'name'
   case fields
   otherwise
    error('Invalid option structure array');
  end
  switch 'type'
   case fields
   otherwise
    [opt.type] = deal([]);
  end
  switch 'index'
   case fields
   otherwise
    [opt.index] = deal([]);
  end
  switch 'predicate'
   case fields
   otherwise
    [opt.predicate] = deal([]);
  end
  switch 'choices'
   case fields
   otherwise
    [opt.choices] = deal([]);
  end
  switch 'values'
   case fields
   otherwise
    [opt.values] = deal([]);
  end
  switch 'value'
   case fields
   otherwise
    [opt.value] = deal([]);
  end
  switch 'assign'
   case fields
   otherwise
    [opt.assign] = deal([]);
  end

  % Default is to take an argument.
  [opt.arg] = deal(1);

  % Type description for error messages.
  [opt.str] = deal('');

  % Predicate function template.
  [opt.p] = deal([]);

  % Check field values.
  for ind = 1:numel(opt)
    if ~ getopt_option_p(opt(ind).name)
      error(sprintf('Option %d has an invalid option name', ind));
    end
    if isempty(opt(ind).type)
      opt(ind).type = 'flag';
    end
    switch opt(ind).type
     case 'alias'
      opt(ind).arg = [];
      if isempty(opt(ind).index)
        error(sprintf('Option %d requires an option index', ind));
      end
     case {'flag', 'counter', 'toggle'}
      opt(ind).arg = 0;
      if isempty(opt(ind).value)
        opt(ind).value = 0;
      elseif ~ isnumeric(opt(ind).value)
        error(sprintf('Option %d requires a numeric option value', ind));
      end
     otherwise
      if isempty(opt(ind).predicate)
        switch opt(ind).type
         case 'color'
          opt(ind).str = 'Matlab color specification';
          opt(ind).p = {2, 'getopt_color_p', []};
         case 'cursor'
          opt(ind).str = 'Matlab cursor name';
          opt(ind).p = {2, 'getopt_cursor_p', []};
         case 'line'
          opt(ind).str = 'Matlab line style';
          opt(ind).p = {2, 'getopt_line_p', []};
         case 'marker'
          opt(ind).str = 'Matlab marker symbol';
          opt(ind).p = {2, 'getopt_marker_p', []};
         case 'symbol'
          opt(ind).str = 'Matlab symbol name';
          opt(ind).p = {1, 'isvarname', []};
         case 'number'
          opt(ind).str = 'number';
          opt(ind).p = {1, 'getopt_number_p', []};
         case 'positive_number'
          opt(ind).str = 'positive number';
          opt(ind).p = {1, 'getopt_positive_number_p', []};
         case 'negative_number'
          opt(ind).str = 'negative number';
          opt(ind).p = {1, 'getopt_negative_number_p', []};
         case 'non_positive_number'
          opt(ind).str = 'non-positive number';
          opt(ind).p = {1, 'getopt_non_positive_number_p', []};
         case 'non_negative_number'
          opt(ind).str = 'non-negative number';
          opt(ind).p = {1, 'getopt_non_negative_number_p', []};
         case 'non_zero_number'
          opt(ind).str = 'non-zero number';
          opt(ind).p = {1, 'getopt_non_zero_number_p', []};
         case 'integer'
          opt(ind).str = 'integral number';
          opt(ind).p = {1, 'getopt_integer_p', []};
         case 'positive_integer'
          opt(ind).str = 'positive integral number';
          opt(ind).p = {1, 'getopt_positive_integer_p', []};
         case 'negative_integer'
          opt(ind).str = 'negative integral number';
          opt(ind).p = {1, 'getopt_negative_integer_p', []};
         case 'non_positive_integer'
          opt(ind).str = 'non-positive integral number';
          opt(ind).p = {1, 'getopt_non_positive_integer_p', []};
         case 'non_negative_integer'
          opt(ind).str = 'non-negative integral number';
          opt(ind).p = {1, 'getopt_non_negative_integer_p', []};
         case 'non_zero_integer'
          opt(ind).str = 'non-zero integral number';
          opt(ind).p = {1, 'getopt_non_zero_integer_p', []};
         case 'any'
          opt(ind).str = 'argument of any type';
          opt(ind).p = {1, 'getopt_any_p', []};
         otherwise
          opt(ind).str = sprintf('argument of class ''%s''', opt(ind).type);
          opt(ind).p = {1, 'isa', [], opt(ind).type};
        end
      else
        opt(ind).str = opt(ind).type;
        switch class(opt(ind).predicate)
         case {'char', 'function_handle'}
          opt(ind).p = {[], opt(ind).predicate, []};
         case 'cell'
          opt(ind).p = {[], opt(ind).predicate{1}, [], opt(ind).predicate{2:end}};
         otherwise
          error(sprintf('Option %d has an invalid option type predicate function', ind));
        end
        % Determine number of return values.
        % This doesn't work for all functions.
        try
          switch nargout(opt(ind).p{2})
           case 1
            opt(ind).p{1} = 1;
           case 2
            opt(ind).p{1} = 2;
           otherwise
            error('Predicate function of option %d has the wrong number of return values', ind);
          end
        end
      end
    end
    if ~ isempty(opt(ind).index)
      if ~ getopt_positive_integer_p(opt(ind).index) | opt(ind).index > numel(opt)
        error(sprintf('Option %d has an invalid option index', ind));
      end
    end
    if ~ isempty(opt(ind).choices) & ~ isempty(opt(ind).values)
      if numel(opt(ind).choices) ~= numel(opt(ind).values)
        error(sprintf('Option %d has a non-matching number of choice/value pairs', ind));
      end
    end
    if ~ isempty(opt(ind).assign)
      if isvarname(opt(ind).assign)
        0;
      elseif opt(ind).assign
        var = opt(ind).name(2:end);
        var(var ~= '_' ...
            & (var < '0' | var > '9') ...
            & (var < 'A' | var > 'Z') ...
            & (var < 'a' | var > 'z')) = '_';
        if ~ isvarname(var)
          error(sprintf('Can not derive variable name for option %d', ind));
        end
        opt(ind).assign = var;
      end
    end
  end

  % Check alias loops.
  for ind = 1:numel(opt)
    switch opt(ind).type
     case 'alias'
      next = ind;
      while strcmp(opt(next).type, 'alias')
        next = opt(next).index;
        if next == ind
          error(sprintf('Option %d refers to itself', ind));
        end
      end
    end
  end

  % Known option names.
  options = {opt.name};

  % Index of the next element of the argument list to be processed.
  optind = 1;

  % While scanning the argument list, options and option arguments
  % are removed from the list.  At the end, the list contains all
  % the remaining arguments.
  while optind <= numel(varargin)

    % Next element of the argument list.
    optopt = varargin{optind};

    % The special argument '--' terminates parsing in all cases.
    if ischar(optopt) & strcmp(optopt, '--')
      varargin(optind) = [];
      break;
    end

    % Skip over non-option arguments.
    if ~ getopt_option_p(optopt)
      optind = optind + 1;
      continue;
    end

    % Remove option from argument list.
    varargin(optind) = [];

    % Check for embedded option argument.
    has_arg = strfind(optopt, '=');
    if has_arg
      start = has_arg(1);
      % Extract option argument (may be empty, take care).
      optarg = optopt(start:end);
      optarg(1) = [];
      % Evaluate simple (constant) arguments.
      optarg = getopt_eval(optarg);
      % Delete option argument from option.
      optopt(start:end) = [];
    end

    % Only process known options.
    switch optopt
     case options
     otherwise
      error(sprintf('Unrecognized option ''%s''', optopt));
    end

    % Determine option structure index (can't fail).
    for ind = 1:numel(opt)
      if strcmp(optopt, opt(ind).name)
        break;
      end
    end

    % Resolve aliases.
    switch opt(ind).type
     case 'alias'
      optval = opt(ind).value;
      while strcmp(opt(ind).type, 'alias')
        ind = opt(ind).index;
        if isempty(optval)
          optval = opt(ind).value;
        end
      end
     otherwise
      optval = [];
    end

    % Update the option value.
    if opt(ind).arg == 0
      % Option has no argument.
      if has_arg
        error(sprintf('Option ''%s'' does not take an argument', optopt));
      end
      if isempty(optval)
        switch opt(ind).type
         case 'flag'
          optval = 1;
         case 'counter'
          optval = opt(ind).value + 1;
         case 'toggle'
          optval = ~ opt(ind).value;
         otherwise
          error(sprintf('Option %d has an unknown option type', ind));
        end
      end
    else
      % Option has an argument.  The argument is either provided implicitly
      % by the option value of an alias option or explicitly by an option
      % argument.
      if has_arg
        if ~ isempty(optval)
          % Option argument supplied via '-opt=arg', but option value
          % is provided through an alias option.
          error(sprintf('Option ''%s'' does not take an argument', optopt));
        end
      else
        if ~ isempty(optval)
          % Option value is provided through an alias option.
          optarg = optval;
        else
          if numel(varargin) < optind
            error(sprintf('Option ''%s'' requires an argument', optopt));
          end
          % Get option argument from the argument list.
          optarg = varargin{optind};
          varargin(optind) = [];
          if getopt_option_p(optarg)
            error(sprintf('Option ''%s'' requires an argument', optopt));
          end
        end
      end
      if isempty(opt(ind).p)
        optval = optarg;
      else
        % Call predicate function.
        opt(ind).p{3} = optarg;
        switch opt(ind).p{1}
         case 1
          true = feval(opt(ind).p{2:end});
          optval = optarg;
         case 2
          [true, optval] = feval(opt(ind).p{2:end});
         otherwise
          try
            [true, optval] = feval(opt(ind).p{2:end});
          catch
            true = feval(opt(ind).p{2:end});
            optval = optarg;
          end
        end
        if ~ true
          switch lower(opt(ind).type(1))
           case {'a', 'e', 'i', 'o', 'u'}
            str = ['an ',  opt(ind).str];
           otherwise
            str = ['a ',  opt(ind).str];
          end
          error(sprintf('Option ''%s'' requires %s', optopt, str));
        end
      end
      if ~ isempty(opt(ind).choices)
        tem = 0;
        for i = 1:numel(opt(ind).choices)
          if iscell(opt(ind).choices)
            elem = opt(ind).choices{i};
          else
            elem = opt(ind).choices(i);
          end
          if isequal(optarg, elem)
            tem = i;
            break;
          end
        end
        if tem == 0
          error(sprintf('Invalid argument for option ''%s''', optopt));
        end
        if ~ isempty(opt(ind).values)
          if iscell(opt(ind).values)
            optval = opt(ind).values{i};
          else
            optval = opt(ind).values(i);
          end
        end
      end
    end

    % Assign option value.
    opt(ind).value = optval;
  end

  % Remove undocumented option structure fields.
  opt = rmfield(opt, {'arg', 'str', 'p'});

  % Remaining arguments.
  rest = varargin;

  % Propagate option values to the caller's workspace.
  for ind = 1:numel(opt)
    if opt(ind).assign
      assignin('caller', opt(ind).assign, opt(ind).value);
    end
  end

% Attempt to evaluate argument ARG.
function val = getopt_eval(arg)

  if exist(arg)
    val = arg;
  else
    val = eval(arg, 'arg');
  end

% Return non-zero if ARG looks like an option.
function p = getopt_option_p(arg)

  if ischar(arg) & numel(arg) > 1
    p = (arg(1) == '-');
  else
    p = 0;
  end

% Return non-zero if ARG is a valid color.
function [p, val] = getopt_color_p(arg)

  val = [];

  switch class(arg)
   case 'char'
    switch arg
     case {'y', 'yellow', ...
           'm', 'magenta', ...
           'c', 'cyan', ...
           'r', 'red', ...
           'g', 'green', ...
           'b', 'blue', ...
           'w', 'white', ...
           'k', 'black'};
      val = arg;
    end
   case 'double'
    switch numel(arg)
     case 1
      if arg >= 0 & arg <= 1
        val = repmat(arg, 1, 3);
      end
     case 3
      if all(arg >= 0 & arg <= 1)
        val = arg;
      end
    end
  end

  p = ~ isempty(val);

% Return non-zero if ARG is a valid cursor.
function [p, val] = getopt_cursor_p(arg)

  val = [];

  switch class(arg)
   case 'char'
    switch arg
     case {'crosshair', 'arrow', 'watch', 'topl', 'topr', 'botl', 'botr', ...
           'circle', 'cross', 'fleur', 'left', 'right', 'top', 'bottom', ...
           'fullcrosshair', 'ibeam'};
      val = arg;
    end
  end

  p = ~ isempty(val);

% Return non-zero if ARG is a valid line style.
function [p, val] = getopt_line_p(arg)

  val = [];

  if isempty(arg)
    arg = 'none';
  end

  switch class(arg)
   case 'char'
    switch arg
     case {'none', ...
           '-', '--', ':', '-.'}
      val = arg;
     case 'solid'
      val = '-';
     case 'dashed'
      val = '--';
     case 'dotted'
      val = ':';
     case 'dashdot'
      val = '-.';
    end
  end

  p = ~ isempty(val);

% Return non-zero if ARG is a valid marker.
function [p, val] = getopt_marker_p(arg)

  val = [];

  if isempty(arg)
    arg = 'none';
  end

  switch class(arg)
   case 'char'
    switch arg
     case {'none', ...
           '+', 'o', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'}
      val = arg;
     case 'plus'
      val = '+';
     case 'circle'
      val = 'o';
     case {'asterisk', 'star'}
      val = '*';
     case 'point'
      val = '.';
     case {'cross', 'x-mark'}
      val = 'x';
     case 'square'
      val = 's';
     case 'diamond'
      val = 'd';
     case 'up'
      val = '^';
     case 'down'
      val = 'v';
     case 'left'
      val = '<';
     case 'right'
      val = '>';
     case 'pentagram'
      val = 'p';
     case 'hexagram'
      val = 'h';
    end
  end

  p = ~ isempty(val);

% Return non-zero if ARG is a number.
function p = getopt_number_p(arg)

  if isnumeric(arg) & numel(arg) == 1 & imag(arg) == 0
    p = 1;
  else
    p = 0;
  end

% Return non-zero if ARG is a positive number.
function p = getopt_positive_number_p(arg)

  if getopt_number_p(arg) & arg > 0
    p = 1;
  else
    p = 0;
  end

% Return non-zero if ARG is a negative number.
function p = getopt_negative_number_p(arg)

  if getopt_number_p(arg) & arg < 0
    p = 1;
  else
    p = 0;
  end

% Return non-zero if ARG is a non-positive number.
function p = getopt_non_positive_number_p(arg)

  if getopt_number_p(arg) & arg <= 0
    p = 1;
  else
    p = 0;
  end

% Return non-zero if ARG is a non-negative number.
function p = getopt_non_negative_number_p(arg)

  if getopt_number_p(arg) & arg >= 0
    p = 1;
  else
    p = 0;
  end

% Return non-zero if ARG is a non-zero number.
function p = getopt_non_zero_number_p(arg)

  if getopt_number_p(arg) & arg ~= 0
    p = 1;
  else
    p = 0;
  end

% Return non-zero if ARG is an integral number.
function p = getopt_integer_p(arg)

  if getopt_number_p(arg) & rem(arg, 1) == 0
    p = 1;
  else
    p = 0;
  end

% Return non-zero if ARG is a positive integral number.
function p = getopt_positive_integer_p(arg)

  if getopt_integer_p(arg) & arg > 0
    p = 1;
  else
    p = 0;
  end

% Return non-zero if ARG is a negative integral number.
function p = getopt_negative_integer_p(arg)

  if getopt_integer_p(arg) & arg < 0
    p = 1;
  else
    p = 0;
  end

% Return non-zero if ARG is a non-positive integral number.
function p = getopt_non_positive_integer_p(arg)

  if getopt_integer_p(arg) & arg <= 0
    p = 1;
  else
    p = 0;
  end

% Return non-zero if ARG is a non-negative integral number.
function p = getopt_non_negative_integer_p(arg)

  if getopt_integer_p(arg) & arg >= 0
    p = 1;
  else
    p = 0;
  end

% Return non-zero if ARG is a non-zero integral number.
function p = getopt_non_zero_integer_p(arg)

  if getopt_integer_p(arg) & arg ~= 0
    p = 1;
  else
    p = 0;
  end

% Return non-zero if ARG is an argument of any type.
function p = getopt_any_p(arg)

  p = 1;

% local variables:
% time-stamp-line-limit: 200
% time-stamp-format: "%:y-%02m-%02d %02H:%02M:%02S %Z"
% page-delimiter: "^function\\>.*\n?"
% end:

%% getopt.m ends here
