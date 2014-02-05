function st = cfu_get_psf_metrics(varargin)
  %function st = cfu_get_psf_metrics(options)
  %
  % Select a region using the mouse. A struct containing min, max,
  % mean, cystic resolution, and FWHM is returned
  %
  % options:
  %       'mode'      'circle' or 'box'                          (default: 'box')
  %       'coord'     display coordinates                            (default: 1)
  %       'log'       logarithmic data, 20*log(ratio)                (default: 1)
  %       'nover'     oversampling used for FWHM interpolation       (default: 4)
  %       'script'    enable scripting, rectangle option must be set (default: 0)
  %       'rect'      rectangle used for inspection, when script is enabled.
  %                   Definition: [UL-x UL-y Span-x Span-y], UL = upper left.
  %       'interp'    interpolation method used. Accepts the same methods as 
  %                   interp1.                                (default: 'spline')
  %       'show'      show ROI used
  %       'fh'        figure handle to grab data from.
  %       'ctr'       Get clutter energy total energy ratio.      (default: true)
  %       'clutter'   Get clutter as function of radius.         (default: false)
  %       'plot'      Plots the CTR curve. Only has an effect when 'ctr' is also true.
  %
  %
  
  if (nargin == 1) && strcmp(varargin{1}, 'test')
    st = inspect_fwhm_test;
    return;
  end
  
  % Defaults
  opt.mode              = 'box';%circle'; % Check for this
  opt.coord             = 1;
  opt.log               = 1;
  opt.script            = 0;
  opt.rect              = [];
  opt.nover             = 4;
  opt.interp            = 'spline';
  opt.show              = 0;
  opt.fh                = [];
  opt.ctr               = true;
  opt.clutter           = false;
  opt.ctr_limit         = 20; % dB
  opt.plot              = false;
  opt.debug             = false;
  opt = va_arg(opt, varargin);
  
  if ~opt.script
    fprintf('Select region using the mouse\n')
    if opt.coord
      rect = rubberband('-return=rectangle',sprintf('-mode=%s',opt.mode),...
                        '-anim=xor','-coord');
    else
      rect = rubberband('-return=rectangle',sprintf('-mode=%s',opt.mode),...
                        '-anim=xor');
    end
  else
    rect = opt.rect;
  end
  
  if isempty(opt.fh)
    hFigure = get(0, 'CurrentFigure');
  else
    hFigure = opt.fh;
  end
  
  if isempty(hFigure)
    error('Can not determine current figure');
  end

  % Original 'NextPlot', 'CurrentAxes'
  opt.CurrentAxis = get(hFigure,'CurrentAxes');
  
  % Find axis for non-Colorbar image
  tem = get(hFigure, 'Children');
  for i=1:length(tem)
    if ~strcmp(get(tem(i),'Tag'),'Colorbar')
      hAxis = tem(i);
    end
  end
  
  st = inspect_fwhm_do(rect,hAxis,opt);
  set(hFigure,'CurrentAxes',opt.CurrentAxis);
  
  
  % Plot
  if opt.ctr && opt.plot
      figure;
      plot(st.radius, st.ct, 'LineWidth',2)
      xlabel('Cyst radius [mm]')
      ylabel('Relative Intensity [dB]')
      
      if st.radius20dB > -Inf,
          hold on;
          plot(st.radius20dB, -20, 'ro', 'LineWidth',2, 'MarkerSize',10)
          hold off;
      end
  end
end






function st = inspect_fwhm_do(rect,hAxis,opt)
  st.rect = rect;
  hDClist = get(hAxis,'Children');
  
  hDC = hDClist(1);
  
  for i=1:length(hDClist)
    if strcmp(get(hDClist(i),'Type'),'image')
      hDC = hDClist(i);
      break;
    end
  end
  
  if ~strcmp(get(hDC,'Type'),'image')
    error('Figure must contain an image')
  end
  

    cdata = get(hDC,'CData');
  % remove -Inf
  if min(cdata(:)) == -Inf,  
      real_data = cdata(cdata>-Inf);
      real_min = min(real_data(:));
      cdata(cdata==-Inf) = real_min;
  end
  
  
  % Size in pixels
  [ny nx] = size(cdata);

  
  % Coordinates on axes
  xdata = get(hDC,'XData');
  ydata = get(hDC,'YData');
  
  xdata = linspace(xdata(1), xdata(end), nx);
  ydata = linspace(ydata(1), ydata(end), ny);
  
  % Range on axes
  xrange = xdata(end)-xdata(1);
  yrange = ydata(end)-ydata(1);
  
  
  % convert values from rect var
  x_min_in  = rect(1);
  y_min_in  = rect(2);
  x_span_in = rect(3);
  y_span_in = rect(4);


  %find closest data points to rect variables
  [dummy roi_x_min_idx] = min( (xdata - x_min_in).^2);
  [dummy roi_y_min_idx] = min( (ydata - y_min_in).^2);
  [dummy roi_x_max_idx] = min( (xdata - (x_min_in+x_span_in)).^2);
  [dummy roi_y_max_idx] = min( (ydata - (y_min_in+y_span_in)).^2);
  
  

  
  % Rectangular contour surrounding region - pixels
  x = nx*(rect(1) - xdata(1))/xrange;
  y = ny*(rect(2) - ydata(1))/yrange;
  
  dx = nx*rect(3)/xrange;
  dy = ny*rect(4)/yrange;
  
  if (strcmp(opt.mode,'circle'))
    hEllipse = imellipse(hAxis,rect);
    api = iptgetapi(hEllipse);
    api.setResizable(false);
    fcn = makeConstrainToRectFcn('imellipse',rect(1)+[0 rect(3)], rect(2)+[0 ...
                        rect(4)]);
    setPositionConstraintFcn(hEllipse,fcn);
    
    tmp = getVertices(hEllipse);
    r = ny*(tmp(:,2)-ydata(1))./yrange;
    c = nx*(tmp(:,1)-xdata(1))./xrange;
    
    if ~opt.show
      delete(hEllipse);
    end
  elseif strcmp(opt.mode,'box')

    % Contour surrounding region (UL, LL, LR, UR)
    r = [floor(y) ceil(y+dy)  ceil(y+dy) floor(y)];
    c = [floor(x) floor(x)   ceil(x+dx)  ceil(x+dx)];
    
    if opt.show
      hRect = imrect(hAxis,rect);
      api = iptgetapi(hRect);
      api.setResizable(false);
      fcn = makeConstrainToRectFcn('imrect',rect(1)+[0 rect(3)], rect(2)+[0 ...
                          rect(4)]);
      setPositionConstraintFcn(hRect,fcn);
    end
  end
  
  
  % get chosen area (ROI) 
  if(roi_x_max_idx > roi_x_min_idx)
    area = cdata(roi_y_min_idx:roi_y_max_idx , roi_x_min_idx:roi_x_max_idx);
  else
    area = cdata(roi_y_min_idx:roi_y_max_idx , roi_x_max_idx:roi_x_min_idx);
  end
  
  % Get FWHM
  if opt.log
      detail_resolution = fwhm_area (area, xdata, ydata, opt.debug);
  else
      detail_resolution = fwhm_area (20*log10(area), xdata, ydata, opt.debug);
  end

  
  st.fwhm_x     = detail_resolution.fwhm_x;
  st.fwhm_y     = detail_resolution.fwhm_y;
  st.fwhm_max   = detail_resolution.fwhm_max_val;
  st.detail_res = detail_resolution;
  
  
  
  
  
  % Mask for extracting region
  mask = roipoly(cdata,c,r);
  
 
  % Data inside mask
  roi =  cdata(mask(:));
  [st.max imax] = max(roi(:));

  % Global indices - TODO: Fix xy issue
  j = 1:size(cdata,1);
  i = 1:size(cdata,2);
  
  js = ndgrid(j,i);
  is = meshgrid(i,j);

  % Global indices of masked pixels (inside selected region)
  is = is(mask(:));
  js = js(mask(:));
  
  % Position of maximum (one-dimensional)
  subixmax = is(imax);
  subiymax = js(imax);
  
  % Position of maximum
  st.x_coord = subixmax*xrange/nx + xdata(1);
  st.y_coord = subiymax*yrange/ny + ydata(1);

  
  
  if opt.log
    % Not necessary to subtract min(roi(:))
    cst.total_energy = sum(10.^(roi./10));
  else
    cst.total_energy = sum(roi.^2);
  end
  
  % Adjust nover to have equal axial and lateral resolution
  if (dx==max(dx,dy))
    nover_x = opt.nover;
    nover_y = ceil(opt.nover*dy/dx);
  else
    nover_y = opt.nover;
    nover_x = ceil(opt.nover*dx/dy);
  end
  % Hack to fix nover_ is sometimes negative
  nover_x = abs(nover_x);
  nover_y = abs(nover_y);
  
  cst.total_energy = cst.total_energy * nover_x * nover_y;  
  %cst.max_radius = min(rect(3),rect(4))/2;
  cst.max_radius = max(rect(3),rect(4))/2;

  cst.cdata = imresize(cdata,[nover_y*ny nover_x*nx]);
  
  if opt.log
      % revert to linear scale
      cst.cdata = 10.^(cst.cdata./20);
  end
  
  % Power
  cst.cdata = cst.cdata.^2;
  
  [cst.ny cst.nx]  = size(cst.cdata);
  
  i = 1:cst.nx;
  j = 1:cst.ny;
  cst.js = ndgrid(j,i);
  cst.is = meshgrid(i,j);

  cst.imax = (st.x_coord - xdata(1))*cst.nx/xrange;
  cst.jmax = (st.y_coord - ydata(1))*cst.ny/yrange;
  cst.xrange = xrange;
  cst.yrange = yrange;
  cst.ctr_limit = opt.ctr_limit;
  

  
  % TODO: Less oversampling  
  if (opt.ctr)
    st.radius = (1:100)/100*cst.max_radius;
  
    % TEST: Works with powers
    mask = imresize(mask,[nover_y*ny nover_x*nx]);
    roi =  cst.cdata(mask(:));
    cst.total_energy = sum(roi);
    
    if opt.debug, figure;imagesc(20*log10(mask.*cst.cdata)); end;

    % Mask points inside ellipsis
    c1 = (cst.is - cst.imax).^2;
    c2 = (cst.js - cst.jmax).^2;
    c3 = cst.max_radius*cst.nx/cst.xrange;
    c4 = cst.max_radius*cst.ny/cst.yrange;
    for i=1:100
        if opt.debug, fprintf('%i ', i); end
        beta = i/100;
        imask = (c1)/ (beta*c3)^2 + (c2)/ (beta*c4)^2 < 1;
        erg_roi =  sum(cst.cdata(imask(:)));
        st.ct(i) = 10*log10(1-abs(erg_roi/cst.total_energy));
    end
    if opt.debug, fprintf('\n'); end
    
    % make sure CT curve contains no complex values
    if ~isreal(st.ct)
        complex_sample = [];
        for idx=1:length(st.ct)
            if ~isreal(st.ct(idx))
                complex_sample(end+1) = idx;
            end
        end
        st.ct(complex_sample)     = [];
        st.radius(complex_sample) = [];

        warning_str = ['CT contained complex values at idx: ' mat2str(complex_sample)];
        %warning(warning_str)
        fprintf(2, '%s\n', warning_str);
    end

    
    % make sure CT curve contains no -Inf values
    if sum(st.ct == -Inf) > 0
        inf_sample = [];
        for idx=1:length(st.ct)
            if st.ct(idx) == -Inf
                inf_sample(end+1) = idx;
            end
        end
        st.ct(inf_sample)     = [];
        st.radius(inf_sample) = [];

        warning_str = ['CT contained -Inf values at idx: ' mat2str(inf_sample)];
        %warning(warning_str)             % print to stdout
        fprintf(2, '%s\n', warning_str); % print to stderr (log)
    end

    
    % make sure CT curve contains no repetitions
    if length(unique(st.ct)) ~= length(st.ct)
        idx_ar = 1:length(st.ct);
        % only keep the unique values
        [dummy uniq_idx] = unique(st.ct);
        %find the missing idx's
        idx_flag = ismember(idx_ar, uniq_idx);
        idx_mis  = idx_ar(~idx_flag);
        
        ct     = st.ct(idx_flag);
        radius = st.radius(idx_flag);
               
        warning_str = ['CT contained repetitive values at idx: ' mat2str(idx_mis)];
        %warning(warning_str)             % print to stdout
        fprintf(2, '%s\n', warning_str); % print to stderr (log)
    else
        ct     = st.ct;
        radius = st.radius;
    end

    
    

    %extrapval = max(st.radius);
    error_val = Inf;
    cystic_res = interp1(ct, radius,  [-6 -12 -20 -40], 'spline', error_val);
    st.radius6dB  = cystic_res(1);
    st.radius12dB = cystic_res(2);
    st.radius20dB = cystic_res(3);
    st.radius40dB = cystic_res(4);
  end
  

  if (opt.clutter)
      bla = imresize(cdata,[nover_y*ny nover_x*nx]);
      %xs = linspace(xdata(1),xdata(end),nx*opt.nover);
      xs = linspace(xdata(1),xdata(end),nx*nover_x);
      xs = meshgrid(xs, bla(:,1));
      %ys = linspace(ydata(1),ydata(2),ny*opt.nover);
      ys = linspace(ydata(1),ydata(2),ny*nover_y);
      ys = ndgrid(ys,bla(1,:));
      thetas = 2*pi*[0:100-1]/100;
      for i=1:100
          xi = (i-1)/100*cst.max_radius*cos(thetas)+st.x_coord;
          yi = (i-1)/100*cst.max_radius*sin(thetas)+st.y_coord;
          st.clutter(i) = mean(interp2(xs,ys,bla,xi,yi)) - st.max;
      end
  end
  
  
end




function [ratio] = inspect_fwhm_contrast_over(beta,cst)
  
  % Mask points inside ellipsis
  imask = ((cst.is - cst.imax).^2) / (beta*cst.max_radius*cst.nx/cst.xrange)^2 + ...
      ((cst.js - cst.jmax).^2) / (beta*cst.max_radius*cst.ny/cst.yrange)^2 < 1;

  % Energy inside ellipsis
  erg_roi =  sum(cst.cdata(imask(:)));

  ratio = abs(erg_roi/cst.total_energy - (1-10.^(-cst.ctr_limit/10)));
end

















function st_out = fwhm_area (area, x_axis, y_axis, dbg_on)
%
%  st_out = fwhm_area (area, x_axis, y_axis [, dbg_on])
%
%  Input:
%     area   = beamformed, scanline converted and log-compressed image
%     x_axis = values of x-axis [mm]
%     y_axis = values of y-axis [mm]
%
%   Returns:
%     struct containing min and max FWHM/HWHM/HWTM
%
% TODO: Implement noise filter when to be used on real measurements
%
% 
%   By Morten F. Rasmussen, 
%         2011-03-15  Init version.
%         2011-03-20  Updated to locate center at center of mass.
%         2013-10-13  Updated to work when the resolution is not the same in both
%                     dimensions. Now also calculates FWHM instead of HWHM, in all directions.
%
    
%% -------------------------------------------
% Vars
ang_step = 128;
[dim_y dim_x] = size(area);
dx = x_axis(2) - x_axis(1);
dy = y_axis(2) - y_axis(1);

if nargin < 4, dbg_on = 0; end



%% -------------------------------------------
% Do the Work!

% Normalise amplitude
area = area - max(area(:));

% get pixel coordinate of maximum in the C-Scan
[y0_i,x0_i] = find(area==0);
if length(x0_i) > 1
    [y0_i,x0_i] = find(area>-6);
    x(1,:) =y0_i; % y-coord
    x(2,:) =x0_i; % x-coord

    cm = [0; 0];
    M=0;
    for idx =1:size(x,2),
        cm = cm + x(:,idx) * area(x(1,idx), x(2,idx));
        M  = M + area(x(1,idx), x(2,idx));
    end
    cm = cm/M;
    y0_i = cm(1);
    x0_i = cm(2);
end
    
x0 = interp1(x_axis,x0_i);
y0 = interp1(y_axis,y0_i);

% $$$ [x_max x_max_i]  = max(area, [], 2);
% $$$ [max_val y0_i]   = max(x_max);
% $$$ x0_i = x_max_i(y0_i);  %x0 in pixel index

if (dbg_on)
    figure(1);
    imagesc(x_axis, y_axis,area, [-60 0]);
    hold on
    plot(x0,y0,'ob');
    hold off
end


%----------------------------
% sampling directions
% We want to make sure we sample in x- and y-direction, in order to be able 
% compare with the standard.
ang_x = [0 pi];
ang_y = [pi/2 -pi/2];
ang   = [((1:ang_step)-1)*2*pi/ang_step ang_x ang_y];

% Straight line sample points. Matrix dim: (angle, dist)
samp_i = 1/20; %sample interval in px
samp_max = sqrt(dim_x^2 + dim_y^2);
dist = 0:samp_i:samp_max;
S_x  = x0_i + cos(ang)'*dist;
S_y  = y0_i + sin(ang)'*dist;

%----------------------------
% Do the actual sampling
extrapval = -Inf; %value of samples from outside domain
method = 'linear';
S = interp2(area, S_x, S_y, method, extrapval);

% We want to search from the out side towards the center
S= fliplr(S);
if (dbg_on)
    x_axis2 = dist(end:-1:1)*dx;
    y_axis2 = ang(end-4:-1:1)*180/pi;
    figure(2); 
    imagesc(x_axis2, y_axis2, S(1:end-4,:), [-60 0])
    hold on;
end

db6  = zeros(ang_step,1);
db20 = zeros(ang_step,1);


for line=1:size(S_x,1),
    % Find first index where value is higher than -6dB
    thr = -6;
    db6_tmp = find((S(line,:)>=thr), 1, 'first');
    % make sure there there also was a value lower than 6dB
    if (isempty(find((S(line,:)<=thr), 1, 'first')) || isempty(db6_tmp))
        % result cannot be used.
        db6_tmp=Inf;
    end
    db6(line) = db6_tmp;
    
    
    thr = -20;
    % Find first index where value is higher than -20dB
    db20_tmp = find((S(line,:)>=thr), 1, 'first');
    % make sure there there also was a value lower than 20dB
    if (isempty(find((S(line,:)<=thr), 1, 'first')) || isempty(db20_tmp))
        % result cannot be used.
        db20_tmp=Inf;
    end
    db20(line) = db20_tmp;

end




% Gather values into standard measures 
max_index = size(S_x,2); % is needed because we sampled from the outside and inwards
db6  = (max_index - db6)*samp_i;
db20 = (max_index - db20)*samp_i;
% convert image pixel index to real units 
dr   = sqrt((sin(ang)*dy).^2 + (cos(ang)*dx).^2);
db6  = db6.*dr';
db20 = db20.*dr';


%% hwhm to fwhm
%indexes of oposing angles.
idx = [1:ang_step/2;  ...
       ang_step/2+1:ang_step];

hwhm = db6(idx);
st_out.fwhm_ar = hwhm(1,:) + hwhm(2,:);
% hwtm to fwtm
hwtm = db20(idx);
st_out.fwtm_ar = hwtm(1,:) + hwtm(2,:);

if (dbg_on >0)
    hold on
    plot(db6(1:end-4),  ang(end-4:-1:1)*180/pi, 'ob');
    plot(db20(1:end-4), ang(end-4:-1:1)*180/pi, 'ok');
    colorbar
    legend('-6 dB', '-20 dB')
    hold off
    xlabel('Angle from maximum [deg]')
    ylabel('Sampling direction [deg]')
end

% TODO: set to Inf if distance >= r

% Statistics on FWHM
[st_out.fwhm_min_val  fwhm_min_dir] =  min(st_out.fwhm_ar(:));
[st_out.fwhm_max_val  fwhm_max_dir] =  max(st_out.fwhm_ar(:));
st_out.fwhm_mean  =  mean(st_out.fwhm_ar(:));
st_out.fwhm_std   =  std(st_out.fwhm_ar(:));
st_out.fwhm_nstd  =  st_out.fwhm_std/st_out.fwhm_mean;

% index to angle [rad]
st_out.fwhm_min_dir  = ang(fwhm_min_dir);
st_out.fwhm_max_dir  = ang(fwhm_max_dir);



%Statistics on dB20
[st_out.fwtm_min_val  fwtm_min_dir]  = min(st_out.fwtm_ar(:));
[st_out.fwtm_max_val  fwtm_max_dir]  = max(st_out.fwtm_ar(:));
st_out.fwtm_mean  = mean(st_out.fwtm_ar(:));
st_out.fwtm_std   = std(st_out.fwtm_ar(:));
st_out.fwtm_nstd  = st_out.fwtm_std/st_out.fwtm_mean;

% index to angle [rad]
st_out.hwtm_min_dir   = ang(fwtm_min_dir);
st_out.hwtm_max_dir   = ang(fwtm_max_dir);


% Classic FWHM
st_out.fwhm_x = db6(end-2) + db6(end-3);
st_out.fwhm_y = db6(end)   + db6(end-1);
% Classic FWTM
st_out.fwtm_x = db20(end-2) + db20(end-3);
st_out.fwtm_y = db20(end)   + db20(end-1);

st_out.ang_ar = ang(1:ang_step/2);

end
















function st = inspect_fwhm_test()

  img_size = 300;
  
  x=-img_size:img_size;
  y=-2*img_size:2*img_size;
  
  xs = meshgrid(x,y);
  ys = ndgrid(y,x);

  x = xs;
  y = ys;
  
  sigmas = 25:-5:5;
  
  z = exp((1./(2*sigmas(1)^2)).*(-x.^2 - y.^2)); 
  z = z + exp((1./(2*sigmas(2)^2)).*(-(x-130).^2 - (y-140).^2)); 
  z = z + exp((1./(2*sigmas(3)^2)).*(-(x+130).^2 - (y-140).^2)); 
  z = z + exp((1./(2*sigmas(4)^2)).*(-(x+130).^2 - (y+140).^2)); 
  z = z + exp((1./(2*sigmas(5)^2)).*(-(x-130).^2 - (y+140).^2));
  
  fh = figure;
  imagesc([img_size -img_size],[img_size -img_size],z);
  colormap(gray)
  title('Two-dimenstional Gaussian (\sigma=5.0)');
  colorbar
  
  disp(sprintf('By selection a region containing a peak, you should obtain\n a FWHM of around 2*sqrt(2*log(2))*sigma = (%f, %f, %f, %f, %f)\n\n',2*sqrt(2*log(2))*sigmas(5),...
               2*sqrt(2*log(2))*sigmas(4),2*sqrt(2*log(2))*sigmas(3),2*sqrt(2*log(2))*sigmas(2),2*sqrt(2*log(2))*sigmas(1)));
  
  st = inspect_fwhm('log',0,'show',1,'mode','box');%,'clutter',true);
  pause(1)
  close(fh)
end
