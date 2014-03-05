% Calculate position and direction of every transducer element.
function [x, y, z, dir] = calc_element_position_LINEAR(myClass_object,TotalLayerThickness,pitch,nr_elements_x,nr_elements_y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   TotalLayerThickness     acLayerThickness    
%   pitch                   pitch
%   nr_elements_x           Number of elements in x direction
%   nr_elements_y           Number of elements in y direction
% OUTPUT
%   x
%   y
%   z
%   dir
% DESCRIPTION
%   Calculate position and direction of every transducer element.
%   - First element at positive x-coordinates.
%   - Last element at negative x-coordinates.
%
%   Based on "SM_GenerateElementPositions", 
%   at: \Engine\SourceCode\RTSC\ScanMan\ScanCalc.c
% VERSION		
%   v1  2009-08-24
%   v2  2010-03-11 - Modified argument list to function with transducer
%                    class
% AUTHOR    Jacob kortbek, Martin Christian Hemmsen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Code as before changing to 2D 
% XfirstElement   = pitch*(N-1)/2;
% fRad            = -TotalLayerThickness;
% i = 0;
% 
% while ( i < N )
%     x(i+1)      = XfirstElement - i*pitch;
%     z(i+1)      = fRad;
%     dir(i+1)    = pi/2;
%     % increment counter
%     i = i+1;
% end

%% Code for 2D transducers
XfirstElement   = pitch*(nr_elements_x-1)/2;
YfirstElement   = pitch*(nr_elements_y-1)/2;

fRad            = -TotalLayerThickness;

x = linspace(-XfirstElement,XfirstElement,nr_elements_x);
y = linspace(-YfirstElement,YfirstElement,nr_elements_y);

x = reshape(repmat(x,1,nr_elements_y),[],1);
y = reshape(repmat(y,nr_elements_x,1),[],1);
z = ones(nr_elements_x*nr_elements_y,1)*-TotalLayerThickness;
dir = ones(nr_elements_x*nr_elements_y,1)*pi/2;














