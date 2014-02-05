% Construct multiple instances of Line class in one go.
%
% If any values are missing, they are duplicated, e.g. if only one origin
% is given and multiply directions, lines will be constructed with
% the same origin and multiple directions
%
% obs = bft3_lines(origins, directions, drs, lengths) 
%
% Create line objects
%
% in:
%    origin        float[# lines,3]
%    direction     float[# lines,3]
%    dr            float[# lines,1]
%    length        float[# lines,1]
%
%  or
%
% obs = bft3_lines(options)
%
% options are string-arguments pairs, where the strings equal the
% titles of the inputs above, e.g.
%
%  origin    = [1:4 ; zeros(2,4)]';
%  direction = [zeros(2,4) ; ones(1,4)]';
%  dr        = ones(4,1);
%  len       = 10*ones(4,1);
%  lines = bft3_lines('origin',origin,'direction',direction,...
%                     'dr',dr,'length',len);
%
% $Id: bft3_lines.m,v 1.14 2011-08-02 18:53:51 jmh Exp $
%

%> @file bft3_lines.m
%> @brief Construct multiple line objects in one go
% ======================================================================
%>
%> @brief Function for constructing multiple lines
%>
%> If any values are missing, they are duplicated, e.g. if only one origin
%> is given and multiply directions, lines will be constructed with
%> the same origin and multiple directions
%>
%> @par Calling:
%> \n ob = bft3_line(origin, direction, dr, length) 
%>
%>  or
%> @par
%> \n ob = bf3_apodization(options)        
%>
%> @par Parameters:
%> <table rules="none">
%>  <tr><td>origin    </td><td>float[\#lines,3] </td></tr>
%>  <tr><td>direction </td><td>float[\#lines,3] </td></tr>
%>  <tr><td>dr        </td><td>float      </td></tr>
%>  <tr><td>length    </td><td>float[\#lines]   </td></tr>
%> </table>
%>
%> $Id: bft3_lines.m,v 1.14 2011-08-02 18:53:51 jmh Exp $
%>
function lines = bft3_lines(varargin)
  st.origin    = [0 0 0];
  st.direction = [0 0 1];
  st.dr        = 1;
  st.length    = 100;
  st.viewport  = [];
  
  if (nargin == 1)
    if (strcmp(varargin{1},'test'))
      lines = bft3_lines_test();
    else
      help(mfilename), error(mfilename)
    end
  elseif (strcmp(class(varargin{1}),'char'))
    st = bft3_va_arg(st,varargin);
    if isempty(st.viewport)
      lines = bft3_lines_do(st.origin, st.direction, st.dr, ...
                            st.length);
    else
      lines = bft3_lines_viewport_do(st.viewport);
    end
  elseif (nargin == 4)
    st.origin    = varargin{1};
    st.direction = varargin{2};
    st.dr        = varargin{3};
    st.length    = varargin{4};

    lines = bft3_lines_do(st.origin, st.direction, st.dr, st.length);
  else
    help(mfilename), error(mfilename)
  end
end

%> @brief Internal function
function lines = bft3_lines_viewport_do(viewport)

  if strcmp(viewport.type,'sector')
    line_pitch = viewport.sector / viewport.n_lines;
    thetas = ...
        ((-(viewport.n_lines-1)/2):((viewport.n_lines-1)/2))*line_pitch;

    line_origins = zeros(3,viewport.n_lines);
    line_origins(1,:) = viewport.radius*sin(thetas) + viewport.start_depth*sin(thetas-viewport.angle);
    line_origins(2,:) = zeros(1,viewport.n_lines);
    line_origins(3,:) = -(viewport.radius)* (1-cos(thetas))+viewport.start_depth*cos(thetas-viewport.angle);
    line_origins = line_origins';
    
    line_directions = [sin(thetas-viewport.angle)'  0*sin(thetas-viewport.angle)'  cos(thetas-viewport.angle)'];

    lines = bft3_lines(line_origins,...
                       line_directions,...
                       viewport.dr, ...
                       viewport.line_length);
  else
    line_origins = zeros(3,viewport.n_lines);
    line_origins(1,:) = ((0:viewport.n_lines-1) - ...
                         (viewport.n_lines-1)/2)*viewport.dx;
    line_origins(2,:) = zeros(1,viewport.n_lines);
    line_origins(3,:) = viewport.start_depth*ones(1,viewport.n_lines);

    line_origins = line_origins';
        
    lines = bft3_lines(line_origins, [0 0 1], viewport.dr, viewport.line_length);
    %    error('Unsupported viewport')
  end

end

%> @brief Internal function
function lines = bft3_lines_do(origins, directions, drs, lengths)
  if ~all([size(origins,2) == 3, size(directions,2) == 3,...
           size(drs,2) == 1, size(lengths,2) == 1])
    help(mfilename), error(mfilename);    
  end
  
  n_lines = max([size(origins,1), size(directions,1),...
                 size(drs,1),size(lengths,1)]);

  % Support common origin, directions, dr or length
  if (size(origins,1) == 1)
    origins = repmat(origins,[n_lines 1]);
  end

  if (size(directions,1) == 1)
    directions = repmat(directions,[n_lines 1]);
  end

  if (size(drs,1) == 1)
    drs = repmat(drs,[n_lines 1]);
  end

  if (size(lengths,1) == 1)
    lengths = repmat(lengths,[n_lines 1]);
  end
  
  if ~all([size(origins,1) == n_lines,...
       size(directions,1) == n_lines,...
       size(drs,1) == n_lines,...
       size(lengths,1) == n_lines])
  else
    for i=1:size(origins,1)
      lines(i) = bft3_line(origins(i,:),directions(i,:),drs(i,:), ...
                           lengths(i,:));
    end
  end
end

%> @brief Internal function
function lines = bft3_lines_test()
  origins    = [1:4 ; zeros(2,4)]';
  directions = [zeros(2,4) ; ones(1,4)]';
  drs = ones(4,1);
  lengths = 10*ones(4,1);
  lines = bft3_lines(origins,directions,drs,lengths);
end
