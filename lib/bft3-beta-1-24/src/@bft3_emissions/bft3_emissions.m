% $Id: bft3_emissions.m,v 1.18 2012-04-25 13:05:09 jmh Exp $
classdef bft3_emissions < handle

  % TODO: - Make delays_and_apodization take a vector of origins
  %       - Change positions to either inside (sine) or outside
  %         (tan) the circle
  % Read-only properties    
  %       - Inside or outside should have same theta_pitch
  properties (SetAccess = 'private', GetAccess = 'public')
    name
    radius
    kerf
    pitch
    height
    r_focus
    n_elements
    theta_pitch;
  end

  properties (SetAccess = 'public', GetAccess = 'public')
    c
    n_active_elements
    before_2011_07_14
    f_xmt
    pos
    window
  end
  
  properties (SetAccess = 'private', GetAccess = 'private')
    thetas
  end

  % Data for 8820e
  % n_elements = 192;
  % Field of View = 60; % Degrees (center to center)
  % f0 = 3.5e6;
  % pitch = 0.33e-3;
  % ROC = 60.3e-3;
  % height 13e-3;
  % PZT + Lens + Matching layer(s)
  % d_MLs = 0.304e-3;
  % c_MLs = 2450;
  % d_lens_max = 0.81e-3;
  % c_lens = 1000;
  % r_focus = 65e-3;
  % angular response 22 degrees (equals acceptance angle)
  % theta_pitch = pi*60/180 / 191
  % 2*atan2(pitch/2,radius), 2*asin(pitch/2,radius)
  
  % if Field of View correct + outer radius (tan)
  % tan(pi*60/180 / 191 / 2) * radius * 2 => pitch = 0.330608e-3
  % if Field of View correct + inner radius (tan)
  % sin(pi*60/180 / 191 / 2) * radius * 2 => pitch = 0.330607e-3
  
  % Assume FOV is incorrect
  % guess a pitch, say 0.35e-4 (Field)
  % or assume inner or outer radius
    
  methods
    function obj = bft3_emissions(varargin)
      st.name              = '8820e_corrected';
      st.c                 = 1480;
      st.n_active_elements = 16;
      st.before_2011_07_14 = 0; % Error in acquisition setup
      st.f_xmt             = -0.5;
      st.pitch             = [];
      st.window            = 'Hamming';
      if (nargin == 1)
        if (strcmp(varargin{1},'test'))
          obj = bft3_emissions_test();
          return;
        end
      else
        st = bft3_va_arg(st,varargin);
      end
      
      % obj = jmh_va_arg(obj,varargin); not working for
      % non-public fields
      obj.name = st.name;
      obj.c    = st.c;
      obj.n_active_elements = st.n_active_elements;
      obj.before_2011_07_14 = st.before_2011_07_14;
      obj.window = st.window;
      obj.f_xmt  = st.f_xmt;
      if strncmp(obj.name,'8820e',5)
        if (strcmp(obj.name,'8820e_corrected'))
          obj.radius     = 0.0600; % Wrong but kept for processing
                                   % of old data
        elseif (strcmp(obj.name,'8820e'))
          obj.radius     = 0.0610; % Wrong but kept for processing
                                   % of old data         
        elseif (strncmp(obj.name,'8820e_rcv',9))
          obj.radius     = 60.3e-3 + 0.304e-3 + 0.81e-3;
        elseif (strncmp(obj.name,'8820e_xmt',9) || ...
                strncmp(obj.name,'8820e_elements',14))
          obj.radius     = 60.3e-3;
        else
          obj.radius     = 61e-3;
        end
        
        obj.kerf       = 3.5e-5;
        obj.height     = 13e-3;
        obj.r_focus    = 65/1000;
        obj.pitch      = 3.3e-4;
        obj.n_elements = 192;
        width       = obj.pitch - obj.kerf;

        if ~isempty(st.pitch)
          % Used to get same theta_pitch for inner and outer aperture
          obj.pitch = st.pitch;
        end
        
        if (strcmp(obj.name,'8820e_rcv_inside') || strcmp(obj.name, ...
                                                    '8820e_xmt_inside'))
          obj.theta_pitch = 2*asin(obj.pitch/(2*obj.radius)); % Skip this
        else
          % Like Field II
          obj.theta_pitch = asin(width/obj.radius) + asin(obj.kerf/obj.radius);
        end
        
        obj.thetas      = (((-(obj.n_elements-1)/2):((obj.n_elements-1)/2))*obj.theta_pitch);
        
        obj.pos = ...
            [obj.radius*sin(obj.thetas) ; zeros(1,192) ; -obj.radius*(1-cos(obj.thetas))]';
        
      elseif strcmp(obj.name,'stv4')
        obj.n_elements = 128;
        obj.pitch  = 2.2e-4;
        obj.kerf   = 2.2e-5;
        obj.height = 15e-3; % Wrong
        obj.pos = ...
            [(((-(obj.n_elements-1)/2):((obj.n_elements-1)/2))*obj.pitch)', ...
             zeros(obj.n_elements,2)];
      elseif strcmp(obj.name,'8804')
        obj.n_elements = 192;
        obj.pitch  = 0.208/1000;
        obj.kerf   = 0.035 / 1000;
        obj.height = 4.5/1000;
        obj.r_focus    = 20/1000;
        obj.pos = ...
            [(((-(obj.n_elements-1)/2):((obj.n_elements-1)/2))*obj.pitch)', ...
             zeros(obj.n_elements,2)];
      else
        error('Unsupported aperture');
      end
    end
    
    function pos.set(obj,data) 
      if (size(data,2) == 3)
        pos = data;
      end
    end
    
    function [vs delays apodization] = ...
          delays_and_apodization(obj,varargin)
      st.n_active_elements = obj.n_active_elements;
      st.f_xmt = obj.f_xmt;
      st.origin = [0 0 0];
      st.window = obj.window;
      st.angle  = 0; % Angle with respect to normal on transducer
      st = bft3_va_arg(st,varargin);

      if strncmp(obj.name,'8820e',5)
        % Angle for center of beam
        theta = atan2(st.origin(:,1),...
                                st.origin(:,3)+...
                      repmat(obj.radius,[size(st.origin,1) 1]));

        theta = theta';
        
        if obj.before_2011_07_14
          active_width = obj.pitch*st.n_active_elements;
        else
          active_width = ...
              obj.radius * 2 * sin(obj.theta_pitch* ...
                                   st.n_active_elements/2);
        end
        
        focus_depth = active_width * st.f_xmt;
        
        % Normal to transducer surface
        vs(1,:) = obj.radius*sin(theta) + focus_depth*sin(theta);
        vs(2,:) = 0;
        vs(3,:) = -((obj.radius)* (1-cos(theta)) - ...
                             cos(theta)*focus_depth);

        % Normal to transducer surface (assume normal penetration
        % of lens in transmit) - not better
        % vs(1,:) = obj.radius*sin(theta) + (focus_depth+0.304e-3 + 0.81e-3)*sin(theta);
        % vs(2,:) = 0;
        % vs(3,:) = -((obj.radius)* (1-cos(theta)) - ...
        %             cos(theta)*(focus_depth+0.304e-3 + 0.81e-3));
        
        % Adjust for steering
        vs(1,:) = vs(1,:) - ...
                focus_depth*(sin(theta) - sin(theta-st.angle));
        vs(3,:) = vs(3,:) - ...
                focus_depth*(cos(theta) - cos(theta-st.angle));

        vs = vs';
        
        % Angle for center of apodization
        theta_apo = max(obj.thetas(1)+(st.n_active_elements-1)*obj.theta_pitch/2,...
              theta);
        theta_apo = ...
            min(obj.thetas(end)-(st.n_active_elements-1)*obj.theta_pitch/2,...
              theta_apo);

        theta_apo = theta_apo';
        
        % Center apodization
        center_apo(1,:) = obj.radius*sin(theta_apo);
        center_apo(2,:) = 0;
        center_apo(3,:) = -obj.radius* (1-cos(theta_apo));

        center_apo = center_apo';
        
        % Point on line passing through center_apo and perpendicular
        % to surface
        tmp(1,:) = obj.radius*sin(theta_apo) + sin(theta_apo);
        tmp(2,:) = 0;
        tmp(3,:) = -((obj.radius)* (1-cos(theta_apo)) - ...
                   cos(theta_apo));
        
        % Adjust for steering
        tmp(1,:) = tmp(1,:) - (sin(theta_apo') - ...
                               sin(theta_apo'-st.angle));
        tmp(3,:) = tmp(3,:) - (cos(theta_apo') - ...
                               cos(theta_apo'-st.angle));
        tmp = tmp';
      elseif (strcmp(obj.name,'stv4') || strcmp(obj.name,'8804'))
        active_width = obj.pitch*st.n_active_elements;
        focus_depth = active_width * st.f_xmt;
        
        % Adjust for steering
        vs = st.origin';
        vs(1,:) =  vs(1,:) + sin(-st.angle)*focus_depth;
        vs(3,:) =  vs(3,:) + cos(-st.angle)*focus_depth;
        
        vs = vs';
        
        center_apo = st.origin';
        center_apo(1,:) = max(repmat(obj.pos(1,1)+ ...
                                     (st.n_active_elements-1)*obj.pitch/2, ...
                                     [1 size(st.origin,1)]),center_apo(1,:));
        center_apo(1,:) = min(repmat(obj.pos(end,1)- ...
                                     (st.n_active_elements-1)*obj.pitch/2, ...
                                     [1 size(st.origin,1)]), ...
                              center_apo(1,:));
        center_apo = center_apo';
        
        tmp = st.origin';
        tmp(1,:) = tmp(1,:) + sin(st.angle);
        tmp(3,:) = tmp(3,:) - cos(st.angle);
        
        tmp = tmp';
      end
      
      w_Hamming = @(n)  0.46*cos(pi*n) + 0.54;
      w_Hanning = @(n)  0.50*cos(pi*n) + 0.50;
      w_Rect    = @(n) n<1;
      
      
      
      eval(['w_apodization = w_', st.window,';']);

      for i_emission=1:size(st.origin,1)
      
        % Compute index_norm for apodization based on perpendicular
        % distance from elements to beam 
        index_norm = ...
          sqrt(sum(...
            cross(repmat(center_apo(i_emission,:)-tmp(i_emission,:),[obj.n_elements 1]),...
                  repmat(tmp(i_emission,:),[obj.n_elements 1])-obj.pos)'.^2))';
        
        if obj.before_2011_07_14
          index_norm = index_norm / (st.n_active_elements* ...
                                     obj.pitch/2);
        else
          index_norm = index_norm / (active_width / 2);
        end
        
        
        apodization1 = zeros(1,obj.n_elements);
        apodization1(index_norm < 1) = w_apodization(index_norm(index_norm ...
                                                          < 1));
        % Hack to ensure no more than n_active_elements
        
        inxs = find(apodization1>0);
        while (length(inxs) > st.n_active_elements)
          cinxs = [inxs(1) inxs(end)];
          [dummy ii] = min(apodization1(cinxs));
          apodization1(cinxs(ii)) = 0;
          inxs = find(apodization1>0);
        end
        while (length(inxs) < st.n_active_elements)
          if (inxs(1) > 1)
            apodization1(inxs(1)-1) = eps;
          else
            apodization1(inxs(end)+1) = eps;
          end
          inxs = find(apodization1>0);
        end

        focus  = vs(i_emission,:);
        delays1 = ...
            sqrt(sum(transpose((obj.pos - repmat(focus,[obj.n_elements 1])).^2)))...
            /obj.c;
        
        if (st.f_xmt < 0)
          % Focus behind aperture
          delays1 = max(delays1) - delays1;
          delays1 = max(delays1(find(apodization1>0))) - delays1;
          delays1(find(apodization1==0)) = max(delays1(find(apodization1>0)));
        else
          % Focus in front of aperture
          delays1 = - delays1;
          delays1(find(apodization1==0)) = ...
              min(delays1(find(apodization1>0)));
          delays1 = delays1 - min(delays1);
        end
        apodization(i_emission,:) = apodization1;
        delays(i_emission,:)      = delays1;
      end
    end
  end
end

function obj = bft3_emissions_test()
  obj = bft3_emissions();
end
