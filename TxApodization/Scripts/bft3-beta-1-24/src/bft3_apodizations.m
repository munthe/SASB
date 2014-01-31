% Construct multiple apodizations in one go
%
% If any values are missing, they are duplicated, e.g. if only one reference
% is given and multiply distances, apodizations will be constructed with
% the same reference and multiple distances
%
% objs = bft3_apodizations(aperture, options)
%
% $Id: bft3_apodizations.m,v 1.8 2012-01-19 10:59:54 jmh Exp $

%> @file bft3_apodizations.m
%>
%> @brief Construct multiple apodization objects in one go
%>
%======================================================================
%> @brief Function for constructing multiple apodizations in one go
%>
%> Function for constructing multiple apodizations in one go using the
%> aperture given as the first argument and a number of options given as
%> string-argument pairs, where the strings equal any of the following:
%> ref, distances, values, dynamic, parametric, fixed, window,
%> window_parameter, f, n_active_elements - similar to when constructing
%> single apodization objects.
%>
%> @return array of instances of the bft3_apodization class.
%======================================================================
function apodizations = bft3_apodizations(aperture, varargin)
  % Function for constructing multiple apodizations in one go
  st.ref        = [0 0 1];
  st.distances  = 0;
  st.values     = [];
  st.dynamic    = false;
  st.parametric = true;
  st.fixed   = false;
  st.window  = 'Hamming';
  st.window_parameter = 1.0;
  st.f       = 1;
  st.orientation = [pi/2 0 0];
  st.n_active_elements = uint32(64);

  if (nargin == 1)
    if ischar(aperture)
      if strcmp(aperture,'test')
        apodizations = bft3_apodizations_test();
      else
        help(mfilename), error(mfilename)
      end
    else
      help(mfilename), error(mfilename)
    end
  elseif (nargin == 4)
    st.aperture  = aperture;
    st.ref       = varargin{1};
    st.distances = varargin{2};
    st.values    = varargin{3};
    apodizations = bft3_apodizations_do(st.aperture, st.ref,...
                                        st.distances, st.values,st);
  elseif ischar(varargin{1})
    st = bft3_va_arg(st,varargin);
    st.aperture = aperture;
    if isempty(st.values)
      st.values = ones(size(aperture.pos,1),1);
    end
    apodizations = bft3_apodizations_do(st.aperture, st.ref,...
                                        st.distances, st.values,st);
  else
    help(mfilename), error(mfilename)
  end
end

%> @brief Internal function
function apodizations = bft3_apodizations_do(aperture, ref, distances, values,st)
  if ~all([size(ref,2) == 3, size(values,1) == size(aperture.pos,1)])
    help(mfilename), error(mfilename);    
  end
  
  n_apodizations = max([size(ref,1), size(distances,1),...
                      size(values,3)]);

  % Support common origin, distances, dr or length
  if (size(ref,1) == 1)
    ref = repmat(ref,[n_apodizations 1]);
  end

  if (size(distances,1) == 1)
    distances = repmat(distances,[n_apodizations 1]);
  end
  
  if (size(values,3) == 1)
    values = repmat(values,[1 1 n_apodizations]);
  end

  if ~all([size(ref,1) == n_apodizations,...
       size(distances,1) == n_apodizations,...
       size(values,3) == n_apodizations])
  else
    for i=1:size(ref,1)
      apodizations(i) = bft3_apodization(aperture, ref(i,:), distances(i,:),...
                                         values(:,:,i));
      
      apodizations(i).dynamic           = st.dynamic;
      apodizations(i).parametric        = st.parametric;
      apodizations(i).fixed             = st.fixed;
      apodizations(i).window            = st.window;
      apodizations(i).window_parameter  = st.window_parameter;
      apodizations(i).f                 = st.f;
      apodizations(i).orientation       = st.orientation;
      apodizations(i).n_active_elements = st.n_active_elements;

      
    end
  end
end

function apodizations = bft3_apodizations_test()
  ref    = zeros(2,3);
  distances = [0 1];
  xmt_aperture = bft3_aperture('test');
  values = zeros(size(xmt_aperture.pos,1),length(distances));
  apodizations = bft3_apodizations(xmt_aperture,...
                                   'ref',ref,...
                                   'distances',distances,...
                                   'values',values,...
                                   'dynamic',true);
end
