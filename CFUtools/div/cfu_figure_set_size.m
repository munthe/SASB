function cfu_set_figure_size(FigW, FigH)
%
%  cfu_set_figure_size(FigW, FigH)
% 
% Sets the size of the current figure to: width = FidW, height = FigH.
%
% By MOFI, 2013-09-01, Init version. -trix was found on the internet.
%

% From: http://tex.stackexchange.com/questions/5559/how-to-avoid-large-margins-around-matlab-plot-in-pdf

set(gcf, 'PaperUnits','centimeters','PaperSize',[FigW FigH],...
         'PaperPosition',[0,0,FigW,FigH],'Units','centimeters',...
         'Position',[1,10,FigW,FigH]);
