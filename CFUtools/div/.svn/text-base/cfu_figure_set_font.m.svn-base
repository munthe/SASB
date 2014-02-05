function cfu_set_figure_font(varargin)
%
%  cfu_set_figure_font(key-value pair arguments)
% Arguments accepted:
%  fig_handle -defaults to gcf.
%  font_name  -defaults to 'Times'.
%  font_size  -defaults to 16.
%
% By MOFI, 2013-09-01, Init Version.
%

st.fig_handle = gcf;
st.font_name  ='Times';
st.font_size  = 16;

st = cfu_parse_input_parameters(st, varargin);

st.ax_handle = findall(st.fig_handle,'type','axes');

set(findall(st.fig_handle,'type','text'), 'FontName',st.font_name)
set(findall(st.fig_handle,'type','text'), 'FontSize',st.font_size)

set(st.ax_handle, 'FontName', st.font_name)
set(st.ax_handle, 'FontSize', st.font_size)