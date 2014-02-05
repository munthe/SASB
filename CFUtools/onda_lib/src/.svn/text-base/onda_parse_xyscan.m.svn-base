% This file parses the an Onda SNQ file.
% It assumes the data recorded is a XY scan.
% 
% [data, x_coord,y_coord] = onda_parse_xyscan (file)
%
% By Morten F. Rasmussen, 2011-07-13
%

function [data, x_coord,y_coord] = onda_parse_xyscan (file)
fid = fopen(file);
%goto data matrix
str_seek('[XY Scan Data 0]', fid); 
%get first line = x coords
tline   = fgetl(fid);         
x_coord = parse_line (tline);
% get the rest
i=1;
tline   = fgetl(fid);
while ~isempty(tline)
    line_ar = parse_line (tline);
    data(i,1:length(line_ar)) = line_ar;
    i=i+1;
    tline = fgetl(fid);
end;    
fclose(fid);
% get the y-coord
y_coord = data(:,1);
data = data(:,2:end);
end

function val = parse_line (line_str)
line_str = lower(line_str); %convert 10E-6 to 10e-6
val = str2num(line_str);
end

function str_seek  (str,fid) 
tline=fgetl(fid);
while ~strcmp(str,tline), 
    tline=fgetl(fid);
end;  
end

