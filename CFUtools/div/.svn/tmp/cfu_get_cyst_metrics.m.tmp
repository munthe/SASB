function st = cfu_get_cyst_metrics(varargin)
%function st = cfu_get_cyst_metrics(options)
%
% Select a two regions using the mouse. A struct containing snr, cnr and
% speckle index is returned. The first region is within the cyst, the second region is
% the speckle.
%
% options:
%        'mode'  'circle' or 'box'                                  ('circle')
%         'log'   true or false. Whether or not data is log-compressed. (true)
%       'coord'   display coordinates                                   (true)
%        'show'   show ellipses used                                    (true)
%       'rect1'   Speckle rectangle. 4-element vector.
%       'rect2'   Cyst    rectangle. 4-element vector.
%     'speckle'   Speckle ellipse.   Struct.
%        'cyst'   Cyst ellipse.      Struct.
%
%  'rect1' and 'rect2' can be used both in box-mode and circle-mode. If rect1 and rect2
%  are used in circle-mode, ellipses are made that exactly fill out the two rectangles.
%  The defintion of the rectangle array is: [top-left-x top-left-y, span-x, span-y].
%
%  'speckle' and 'cyst' are only used in circle-mode. They are both a structs
%  containing the following members: 
%     center = [x y], 
%     radius = [x y],  
%
%  If show==1, a handle to both markers (square and ellipse) is output as:
%  h_marker1 and h_marker2. 
%
%  To change the color of a marker call:
%    >> metric = cfu_get_cyst_metrics;
%    >> color = [0 0 0];
%    >> setColor(metric.h_maker1, color);
%
%  To change the line style of a marker call:
%    >> ch = get(metric.h_marker1, 'Children');
%    >> linestyle = '--';
%    >> for idx=1:length(ch), set(ch(idx), 'LineStyle', linestyl); end
%
%  To change the line width cal:
%    >> linewidth = 2;
%    >> for idx=1:length(ch), set(ch(idx), 'LineWidth', linewidth); end
%
%
% TODO: expand the handling / definition of rect1 and rect2 to structs.
%
% Version 1.13 2010/09/14 jmh
% Version 1.14 2013/11/05 MOFI. Added the 'log' parameter. Updated the help
%                               file. Renamed from inspect_cnr to cfu_get_cyst_metrics.
% Version 1.15 2013/11/06 MOFI. Can now also show the marked rectangles.
% Version 1.16 2013/11/06 MOFI. Now also accepts ellipse definition as input.
% Version 1.17 2013/12/18 MOFI. Changed the definition of ellipse 1 and 2. They are now
%                               termed 'speckle' and 'cyst' and are defined by a 1x2
%                               array. Input parameter 'script' no longer needs to be set.
%

if (nargin == 1) && strcmp(varargin{1}, 'test'), inspect_cnr_test, return, end

% Defaults
opt.mode      = 'circle';
opt.script    = [];
opt.rect1     = [];
opt.rect2     = [];
opt.speckle   = [];
opt.cyst      = [];
opt.show      = 1;
opt.log       = 1;

opt = cfu_parse_input_parameters(opt, varargin);

%% Select data
if isempty(opt.script)
    % base parameter 'script' on whether input parameters are set.
    if isempty(opt.speckle) && isempty(opt.speckle) && isempty(opt.rect1) && isempty(opt.rect2) 
        opt.script = false;
    else
        opt.script = true;
    end
end

if ~opt.script 
    % Using the mouse (with "rubberband")
    fprintf('Select first region (speckle) using the mouse - ')
    if strcmp(opt.mode, 'circle')
        fprintf('center then upper right\n')
    else
        fprintf('upper left then lower right\n')
    end
    
    rect1 = rubberband('-return=rectangle',...
                     sprintf('-mode=%s',opt.mode),'-anim=xor','-coord')
    
    fprintf('Select second region (cyst) using the mouse - ')
    if strcmp(opt.mode, 'circle')
        fprintf('center then upper right\n')
    else
        fprintf('upper left then lower right\n')
    end
    rect2 = rubberband('-return=rectangle',...
                       sprintf('-mode=%s',opt.mode),'-anim=xor','-coord')
    
else
    % Scripting
    if strcmp(opt.mode, 'circle') % Circle mode
        if ~isempty(opt.speckle) || ~isempty(opt.cyst)
            % rect1
            rect1_min_x = opt.speckle.center(1) - opt.speckle.radius(1);
            rect1_min_y = opt.speckle.center(2) - opt.speckle.radius(2);
            rect1_span_x = 2*opt.speckle.radius(1);
            rect1_span_y = 2*opt.speckle.radius(2);
            % rect2
            rect2_min_x = opt.cyst.center(1) - opt.cyst.radius(1);
            rect2_min_y = opt.cyst.center(2) - opt.cyst.radius(2);
            rect2_span_x = 2*opt.cyst.radius(1);
            rect2_span_y = 2*opt.cyst.radius(2);
            % gather into array that rubberband accepts
            rect1 = [rect1_min_x, rect1_min_y, rect1_span_x, rect1_span_y];
            rect2 = [rect2_min_x, rect2_min_y, rect2_span_x, rect2_span_y];
            
        else
            warning(['Mode is ''circle'', but the ellipses where not set. Using rect ' ...
                     'difinitions instead.'])
        end
    end
        
    
    if isempty(opt.speckle) % Box mode OR circle mode with "ellipseX" not set
        if length(opt.rect1) ~= 4 || length(opt.rect2) ~= 4
            error('Both rect1 and rect2 must contain 4 elements')
        end
        % The input rectanlges are OK -> use them
        rect1 = opt.rect1;
        rect2 = opt.rect2;
    end
end




hFigure = get(0, 'CurrentFigure');
if isempty(hFigure)
  error('Can not determine current figure');
end

% Find axis for non-Colorbar image
tem = get(hFigure, 'Children');
for i=1:length(tem)
  if ~strcmp(get(tem(i),'Tag'),'Colorbar')
    hAxis = tem(i);
  end
end

st1 = inspect_cnr_do(rect1, hAxis,opt);
st2 = inspect_cnr_do(rect2, hAxis,opt);

<<<<<<< .mine
st.cnr           = (st1.mean - st2.mean) /  sqrt(st1.var + st2.var); % the larger the better
st.speckle_index = sqrt(st1.var)/st1.mean + sqrt(st2.var)/st2.mean;  % the smaller the better
st.snr           = (st1.mean - st2.mean) /  sqrt(st1.mean^2+st2.mean^2); % the larger the better
st.cr            = (st1.mean - st2.mean) / (st1.mean + st2.mean); % the larger the better
=======
st.cnr           = (st1.mean - st2.mean) /  sqrt(st1.var + st2.var);
st.cnr2          = (st1.mean - st2.mean) /  sqrt(st1.var);
st.speckle_index = sqrt(st1.var)/st1.mean + sqrt(st2.var)/st2.mean;
st.snr           = (st1.mean - st2.mean) /  sqrt(st1.mean^2+st2.mean^2);
st.cr            = (st1.mean - st2.mean) / (st1.mean + st2.mean);
>>>>>>> .r82
st.snr0_speckle  = st1.mean / sqrt(st1.var);
st.snr0_cyst     = st2.mean / sqrt(st2.var);

st.rect1 = rect1;
st.rect2 = rect2;


if strcmp(opt.mode, 'circle') % Set output ellipse definition
    if isempty(opt.speckle)
        % ellipse 1
        st.speckle.center(1) = st.rect1(1) + st.rect1(3)/2;
        st.speckle.center(2) = st.rect1(2) + st.rect1(4)/2;
        st.speckle.radius(1) = st.rect1(3)/2;
        st.speckle.radius(2) = st.rect1(4)/2;
        % ellipse 2
        st.cyst.center(1) = st.rect2(1) + st.rect2(3)/2;
        st.cyst.center(2) = st.rect2(2) + st.rect2(4)/2;
        st.cyst.radius(1) = st.rect2(3)/2;
        st.cyst.radius(2) = st.rect2(4)/2;
    else
        st.speckle = opt.speckle;
        st.cyst = opt.cyst;
    end
end

% output marker handles
if opt.show
    st.h_marker1 = st1.h_marker;
    st.h_marker2 = st2.h_marker;    
end

end % end of main function






%% Worker function
function st = inspect_cnr_do(rect,hAxis,opt)

hDClist = get(hAxis,'Children');

hDC = hDClist(1);

for i=1:length(hDClist)
  if strcmp(get(hDClist(i),'Type'),'image')
    hDC = hDClist(i);
    break
  end
end

if ~strcmp(get(hDC,'Type'),'image')
  error('Figure must contain an image')
end

xdata = get(hDC,'XData');
ydata = get(hDC,'YData');
cdata = get(hDC,'CData');

xrange = xdata(end)-xdata(1);
yrange = ydata(end)-ydata(1);
dimc = size(cdata);

if opt.log
    cdata = 10.^(cdata/20);
end

% Rectangular contour surrounding region
x = dimc(2)*(rect(1) - xdata(1))/xrange;
y = dimc(1)*(rect(2) - ydata(1))/yrange;
dx = dimc(2)*rect(3)/xrange;
dy = dimc(1)*rect(4)/yrange;

% image data
r = [floor(y) ceil(y+dy)  ceil(y+dy) floor(y)];
c = [floor(x) floor(x)   ceil(x+dx)  ceil(x+dx)];

if (strcmp(opt.mode,'circle'))
    % Draw elipsis and get its vertices
    hEllipse = imellipse(hAxis,rect);
    api = iptgetapi(hEllipse);
    api.setColor('y');
    api.setResizable(false);
    
    fcn = makeConstrainToRectFcn('imellipse',rect(1)+[0 rect(3)],...
                                 rect(2)+[0 rect(4)]);
    setPositionConstraintFcn(hEllipse,fcn);
    
    tmp = getVertices(hEllipse);
    r = dimc(1)*(tmp(:,2)-ydata(1))./yrange;
    c = dimc(2)*(tmp(:,1)-xdata(1))./xrange;

    f = @(varargin) delete(varargin{1});
    set(hEllipse,'ButtonDownFcn',f)
    set(hEllipse,'HitTestArea','on');

    % delete the ellipse, if show == false
    if ~opt.show
        delete(hEllipse);
    else
        st.h_marker = hEllipse;
    end
elseif (strcmp(opt.mode,'box'))
    % TODO: draw the box

    h_rect = imrect(hAxis,rect);
    api   = iptgetapi(h_rect);
    api.setColor('y');
    api.setResizable(false);
    
    fcn = makeConstrainToRectFcn('imellipse',rect(1)+[0 rect(3)],...
                                 rect(2)+[0 rect(4)]);
    setPositionConstraintFcn(h_rect,fcn);
    
    f = @(varargin) delete(varargin{1});
    set(h_rect,'ButtonDownFcn',f)
    set(h_rect,'HitTestArea','on');

    % delete the box, if show == false
    if ~opt.show
        delete(h_rect);
    else
        st.h_marker = h_rect;
    end
end

% mask
mask = roipoly(cdata,c,r);

%figure;imagesc(mask.*cdata)

% masked data
roi =  cdata(mask(:));

% mean and variance
st.mean = mean(roi(:));
st.var  = var(roi(:));

end









%% Test Function
function inspect_cnr_test()

img_size = 300;

x=-img_size:img_size;

[x y] = ndgrid(x,x);

sigmas = 25:-5:5;

z = exp((1./(2*sigmas(1)^2)).*(-x.^2 - y.^2));
z = z + exp((1./(2*sigmas(2)^2)).*(-(x-130).^2 - (y-130).^2));
z = z + exp((1./(2*sigmas(3)^2)).*(-(x+130).^2 - (y-130).^2)); 
z = z + exp((1./(2*sigmas(4)^2)).*(-(x+130).^2 - (y+130).^2)); 
z = z + exp((1./(2*sigmas(5)^2)).*(-(x-130).^2 - (y+130).^2));

fh = figure;
imagesc([img_size -img_size],[-img_size img_size],z);
colormap(gray)
title('Two-dimenstional Gaussian (\sigma=5.0)');

st = inspect_cnr()
close(fh)

end
