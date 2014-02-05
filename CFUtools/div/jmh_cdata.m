% $Id: jmh_cdata.m,v 1.1 2011/09/10 19:00:41 jmh Exp $
function [cdata xlim ylim] = jmh_cdata(varargin)

  % Defaults
  opt.gca = gca;
  
  opt = cfu_parse_input_parameters(opt,varargin);
  
  hAxis = opt.gca;

  hDClist = get(hAxis,'Children');

  hDC = hDClist(1);

  for i=1:length(hDClist)
    if strcmp(get(hDClist(i),'Type'),'image')
      hDC = hDClist(i);
      break;
    end
  end

  cdata = [];
  xlim  = [];
  ylim  = [];
  
  if strcmp(get(hDC,'CDataMapping'),'scaled')
    xlim = get(hAxis,'XLim');
    ylim = get(hAxis,'YLim');
    xdata = get(hDC,'XData');
    ydata = get(hDC,'YData');
    cdata = get(hDC,'CData');

   
    xrange = xdata(end)-xdata(1);
    yrange = ydata(end)-ydata(1);
    a = round((xlim(1) - xdata(1)) * size(cdata,2) / xrange);
    b = round((ylim(1) - ydata(1)) * size(cdata,1) / yrange);
    a = max(a,1);
    b = max(b,1);
    
    na = floor((xlim(2) - xlim(1)) * size(cdata,2) / xrange);
    nb = floor((ylim(2) - ylim(1)) * size(cdata,1) / yrange);
    nb = min(nb,size(cdata,1)-1);
    na = min(na,size(cdata,2)-1);
    
    cdata = cdata(b:b+nb,a:a+na); % TODO assumes direction of
                                   % y-axis
  end
  
  if min(size(cdata)) == 0
      cdata = get(hDC,'CData');
      warning('This script was failing.. Forced value to ''cdata''.')
  end
end
