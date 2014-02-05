% rubberband -- select a figure or axes region.
%
%      [...] = rubberband(h, ...)
%
% The 'rubberband' command waits until a key or mouse button is pressed
% over a figure and sets the start point at that point.  It then tracks
% the figure's cursor until a key or mouse button is pressed again.  The
% second point marks the end point of the region.
%
% Optional first argument H is the reference axes or figure handle.
% Default is to use the current axes or, if no current axes exist, the
% current figure.  The rest of the argument list are options.  Below is
% a list of all the options, together with their meaning.
%
% -anim ARG
%      Rendering of animated objects.  An argument of 'normal' means to
%      use the double buffer, 'xor' means to draw and erase animated
%      objects with xor bit operations.  The later is faster with the
%      drawback of less accurate colors of animated objects.
% -color ARG
%      Color of animated objects.  Argument ARG is a Matlab color.
%      Default is to use white or black depending on the background
%      color of the reference axes or figure.
% -coord
%      Display cursor coordinates while tracking the figure's cursor.
% -cursor ARG
%      Figure cursor.  Argument ARG is a Matlab cursor name.  Default
%      is to use a cross-hair cursor.  If axes expansion mode is on, the
%      default is to use a full cross-hair cursor.
% -expand ARG
%      Expand axes to figure dimensions.  An argument of 'x' means to
%      expand the x-axis, 'y' means to expand the y-axis, and 'both'
%      means to expand both axes.  If the reference object is a figure,
%      the -expand option has no effect except for selecting the default
%      cursor.
% -line ARG
%      Rubber-band line style.  Argument ARG is a Matlab line style (a
%      string).  Default is to draw a dotted line.  An argument of 'none'
%      means to draw no line.
% -marker ARG
%      Rubber-band marker symbol.  Argument ARG is a Matlab marker symbol
%      (a string).  Default is to draw a point.  An argument of 'none'
%      means to draw no marker.
% -mode ARG
%      Rubber-band drawing mode.  An argument of 'box' means to draw
%      a rectangular box (the default), 'x' means to draw a horizontal
%      line, 'y' means to draw a vertical line, 'line' means to
%      draw a straight line, and 'circle' means to draw an ellipsis
%      inside the box.
% -return ARG
%      Return values.  If argument ARG is 'points', return start point
%      and end point (the default), if ARG is 'vectors', return x-axis
%      values and y-axis values, if ARG is 'rectangle', return rectangle
%      dimensions.
%
% If the -coord option is combined with rectangle dimensions as return
% values, the displayed end point coordinates are the distance between
% the end point and the start point.  Otherwise, cursor coordinates are
% absolute values.
%
% See the 'plot' command, for a list of predefined Matlab colors, line
% styles and marker symbols.
%
% Return value is either a matrix, two vectors, or four scalars.  Size and
% grouping of the return values depends on the -return option according to
% the following table.
%
%  Return values |      points      |      vectors     |     rectangle
% ---------------+------------------+------------------+------------------
%  one matrix    |  [x1 y1; x2 y2]  |  [x1 x2; y1 y2]  |     [x y w h]
%  two vectors   | [x1 y1], [x2 y2] | [x1 x2], [y1 y2] |   [x y], [w h]
%  four scalars  |  x1, y1, x2, y2  |  x1, x2, y1, y2  |    x, y, w, h
%
% For point return values, X1 and Y1 are the start point coordinates,
%  and X2 and Y2 are the end point coordinates.
% For vector return values, X1 and X2 are the x-axis coordinates, and
%  Y1 and Y2 are the corresponding y-axis coordinates.  Coordinates are
%  sorted in ascending order.
% For rectangle return values, X and Y are the coordinates of the corner
%  with the smallest coordinates of the rectangle, and W and H are the 
%  width and height of the rectangle.
%
% When selecting a figure region coordinate values are in figure units.
% Otherwise, coordinate values are in axes data units.

%% rubberband.m --- the preceding comment is the documentation string.

% Copyright (C) 2004 Ralph Schleicher

% Author: Ralph Schleicher <rs@nunatak.allgaeu.org>
% Time-stamp: <2004-07-15 22:49:45 CEST>

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

%% Code:


% Program entry point.
function varargout = rubberband(varargin)

if (nargin == 1) && strcmp(varargin{1}, 'test'), rubberband_test, return, end
  
  
  % Initialize options.
  for ind = 1:intmax('int16')
    switch ind
     case 1
      % Rubber-band mode; 0 means box (the default), 1 means x-axis,
      % 2 means y-axis, and 3 means straight line.
      opt(ind).name = '-mode';
      opt(ind).type = 'char';
      opt(ind).choices = {'box', 'x', 'y', 'line', 'circle'};
      opt(ind).values = 0:4;
      opt(ind).value = 0;
     case 2
      % Line style, default dottet.
      opt(ind).name = '-line';
      opt(ind).type = 'line';
      opt(ind).value = ':';
     case 3
      % Marker symbol, default point.
      opt(ind).name = '-marker';
      opt(ind).type = 'marker';
      opt(ind).value = '.';
     case 4
      % Non-zero means to display cursor coordinates.
      opt(ind).name = '-coord';
     case 5
      % User defined color; default is either black or white.
      opt(ind).name = '-color';
      opt(ind).type = 'color';
     case 6
      % User defined cursor; default is either a cross-hair
      % or full cross-hair cursor.
      opt(ind).name = '-cursor';
      opt(ind).type = 'cursor';
     case 7
      % Axes expansion; 0 means no expansion (the default),
      % 1 means x-axis expansion, 2 means y-axis expansion,
      % and 3 means to expand both axes.
      opt(ind).name = '-expand';
      opt(ind).type = 'char';
      opt(ind).choices = {'x', 'y', 'both'};
      opt(ind).values = 1:3;
      opt(ind).value = 0;
     case 8
      % Return values; 0 means points (the default), 1 means vectors,
      % and 2 means rectangle dimensions.
      opt(ind).name = '-return';
      opt(ind).type = 'char';
      opt(ind).choices = {'points', 'vectors', 'rectangle'};
      opt(ind).values = 0:2;
      opt(ind).value = 0;
     case 9
      % Rendering of animated objects.
      opt(ind).name = '-anim';
      opt(ind).type = 'char';
      opt(ind).choices = {'normal', 'xor'};
      opt(ind).value = 'normal';
     otherwise
      break;
    end
  end

  
  tem = opt;
  
  % Parse options from argument list.
  [opt, arg] = getopt(opt, varargin{:});

  % TODO use global p, f = @(x) p =1; and assign as predicate
  
  % Hack to see if an option is set or just default
  if (opt(1).value == 4)
    % User have selected circle
    if (opt(3).value == '.')
      % Default marker option is active, check if its forced by user
      tem(3).value = 'none';
      [tem, arg] = getopt(tem,varargin{:});
      if (tem(3).value == '.')
        % User have selected default, do nothing
      else
        % Select marker 'none' for circles
        opt(3).value = 'none';
      end        
    end      
  end
  
  % Simplify options structure.
  tem = opt;
  opt = [];

  for ind = 1:numel(tem)
    opt = setfield(opt, tem(ind).name(2:end), tem(ind).value);
  end

  % Assign optional arguments.
  switch numel(arg)
    case 0
      h = get(0, 'CurrentFigure');
      if isempty(h)
        error('Can not determine current figure');
      end
      
      tem = get(h, 'Children');
      for i=1:length(tem)
        if ~strcmp(get(tem(i),'Tag'),'Colorbar')
          h = tem(i);
        end
      end
    case 1
      h = arg{1};
    otherwise
      error(nargchk(0, 1, numel(arg)));
  end
  
  % Get figure handle.
  switch get(h, 'Type')
   case 'figure'
    fig = h;
   case 'axes'
    fig = get(h, 'Parent');
   otherwise
    error('Invalid handle');
  end

  % Save figure properties.
  prop.Pointer = get(fig, 'Pointer');
  prop.Renderer = get(fig, 'Renderer');
  prop.DoubleBuffer = get(fig, 'DoubleBuffer');
  prop.WindowButtonMotionFcn = get(fig, 'WindowButtonMotionFcn');

  % Determine drawing axes dimensions.
  if h == fig
    % Axes position (whole figure).
    pos = [0 0 1 1];

    % Axes limits.
    tem = get(h, 'Position');

    xlim = [0 tem(3)];
    ylim = [0 tem(4)];
    ydir = 'normal';
    xdir = 'normal';
  else
    % Reference axes position.
    units = get(h, 'Units');
    set(h, 'Units', 'normalized');
    pos = get(h, 'Position');
    set(h, 'Units', units);
    ydir = get(h,'YDir');
    xdir = get(h,'XDir');

    % Expand axes limits.
    xlim = get(h, 'XLim');
    if bitand(opt.expand, 1)
      switch get(h, 'XScale')
       case 'log'
        error('Can not expand logarithmic scale');
      end
      xlim = xlim(1) + ([0 1] - pos(1)) ./ pos(3) .* (xlim(2) - xlim(1));
      pos([1 3]) = [0 1];
    end

    ylim = get(h, 'YLim');
    if bitand(opt.expand, 2)
      switch get(h, 'YScale')
       case 'log'
        error('Can not expand logarithmic scale');
      end
      ylim = ylim(1) + ([0 1] - pos(2)) ./ pos(4) .* (ylim(2) - ylim(1));
      pos([2 4]) = [0 1];
    end
  end
  
  % Create drawing axes.
  ax = axes('Parent', fig, ...
            'Visible', 'off', ...
            'Units', 'normalized', ...
            'Position', pos, ...
            'XLim', xlim, ...
            'YLim', ylim, ...
            'YDir', ydir, ...
            'XDir', xdir);

  % Save reference object handle
  % and axes limits.
  if h == fig
    data.ref = fig;
  else
    data.ref = ax;
  end

  data.xlim = xlim;
  data.ylim = ylim;

  % Determine color for animated objects.
  if isempty(opt.color)
    opt.color = get(data.ref, 'Color');
    if ischar(opt.color)
      switch opt.color
       case {'none', 'black', 'k'}
        opt.color = [1 1 1];
       otherwise
        opt.color = [0 0 0];
      end
    elseif all(opt.color < 0.5)
      opt.color = [1 1 1];
    else
      opt.color = [0 0 0];
    end
  end

  % Create rubber-band line object.
  data.line = line('Parent', ax, ...
                   'EraseMode', opt.anim, ...
                   'Color', opt.color, ...
                   'LineStyle', opt.line, ...
                   'Marker', opt.marker, ...
                   'XData', [], ...
                   'YData', []);

  % Optionally create text object to
  % display cursor coordinates.
  if opt.coord
    data.coord = text('Parent', ax, ...
                      'EraseMode', opt.anim, ...
                      'Color', opt.color, ...
                      'HorizontalAlignment', 'left', ...
                      'VerticalAlignment', 'bottom', ...
                      'FontSize', get(ax, 'FontSize') / 1.2);
  else
    data.coord = [];
  end

  % Establish look-up table
  if (opt.mode == 4)
    t = 0:0.01:2*pi;
    cos_t = cos(t);
    sin_t = sin(t);
    data.cos_t = cos_t;
    data.sin_t = sin_t;
  end
  
  % Change figure cursor.
  if isempty(opt.cursor)
    if opt.expand
      opt.cursor = 'fullcrosshair';
    else
      opt.cursor = 'crosshair';
    end
  end

  set(fig, 'Pointer', opt.cursor);

  % Remember if full cross-hair cursor is enabled.
  opt.full = strcmp(opt.cursor, 'fullcrosshair');

  % Adjust figure renderer.
  switch opt.anim
   case 'normal'
    set(fig, 'DoubleBuffer', 'on', ...
             'Renderer', 'painters');
   case 'xor'
    set(fig, 'DoubleBuffer', 'off');
  end

  % Install call-back handler.
  set (fig, 'WindowButtonMotionFcn', {@rubberband_motion, data, opt});

  % Get start point.
  waitforbuttonpress;
  p1 = rubberband_point(data);

  % Draw start point.
  set(data.line, 'XData', p1(1), ...
                 'YData', p1(2), ...
                 'UserData', p1);

  % Get end point.
  waitforbuttonpress;
  p2 = rubberband_point(data);

  % Delete drawing axes.
  delete(ax);

  % Restore figure properties.
  set(fig, 'Pointer', prop.Pointer, ...
           'Renderer', prop.Renderer, ...
           'DoubleBuffer', prop.DoubleBuffer, ...
           'WindowButtonMotionFcn', prop.WindowButtonMotionFcn);

  % Set unneeded values to zero.
  switch opt.mode
    case 1
      p1(2) = 0;
      p2(2) = 0;
    case 2
      p1(1) = 0;
      p2(1) = 0;
    case 4
      % Surrounding rectangle 
      p1 = p1 - (p2-p1);
  end

  % Prepare result.
  switch opt.return
   case 0
    % Return points.
    res = [p1; p2]';
   case 1
    % Return vectors.
    res = sortrows([p1; p2]);
   case 2
    % Return rectangle dimensions.
    res = [min(p1, p2) abs(p1 - p2)]';
  end
  switch nargout
   case {0, 1}
    varargout = {res'};
   case 2
    varargout = {res(1:2), res(3:4)};
   case 4
    varargout = {res(1), res(2), res(3), res(4)};
   otherwise
    error('Wrong number of output arguments');
  end


% Rubber-band call-back handler.
function rubberband_motion(obj, event, data, opt)

  if ~ ishandle(data.line)
    return;
  end

  % Get start point.
  p1 = get(data.line, 'UserData');

  % Get end point.
  p2 = rubberband_point(data);

  % Update rubber-band.
  if ~ isempty(p1)
    if opt.mode == 0
      % Draw rectangular box.
      if opt.full
        x = [p1(1) p1(1) p2(1)];
        y = [p2(2) p1(2) p1(2)];
      else
        x = [p1(1) p1(1) p2(1) p2(1) p1(1)];
        y = [p2(2) p1(2) p1(2) p2(2) p2(2)];
      end
    elseif opt.mode == 4
      xr = abs(p1(1) - p2(1));
      yr = abs(p1(2) - p2(2));
      x = xr*data.cos_t + p1(1);
      y = yr*data.sin_t + p1(2);
    else
      % Draw straight line.
      switch opt.mode
       case 1
        % Horizontal line.
        p2(2) = p1(2);
       case 2
        % Vertical line.
        p2(1) = p1(1);
      end
      x = [p1(1) p2(1)];
      y = [p1(2) p2(2)];
    end
    set(data.line, 'XData', x, ...
                   'YData', y);
  end

  % Update cursor coordinates.
  if ishandle(data.coord)
    if opt.return == 2 & ~ isempty(p1)
      p = p2 - p1;
    else
      p = p2;
    end
    switch opt.mode
     case 1
      str = sprintf(' %.5G', p(1));
     case 2
      str = sprintf(' %.5G', p(2));
     otherwise
      str = sprintf(' %.5G, %.5G', p);
    end
    set(data.coord, 'Position', p2, ...
                    'String', str);
  end


% Return cursor coordinates.
function p = rubberband_point(data)

  % Get coordinates.
  p = get(data.ref, 'CurrentPoint');
  
  switch get(data.ref, 'Type')
    case 'axes'  
      p = p(1, 1:2);
  end

  % Apply limits.
  if p(1) < data.xlim(1)
    p(1) = data.xlim(1);
  end
  if p(1) > data.xlim(2)
    p(1) = data.xlim(2);
  end
  if p(2) < data.ylim(1)
    p(2) = data.ylim(1);
  end
  if p(2) > data.ylim(2)
    p(2) = data.ylim(2);
  end


% Local test procedure.
function rubberband_test

  x = 0:pi/100:2*pi;

  subplot(2, 1, 1);
  plot(x, sin(x));
  axis tight;

  subplot(2, 1, 2);
  plot(x, cos(x));
  axis tight;

  [p1, p2] = rubberband('-mode=x', '-coord', '-expand=y');


%% rubberband.m ends here
