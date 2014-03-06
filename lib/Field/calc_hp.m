%  Procedure for calculating the emitted field.
%
%  Calling:  [hp, start_time] = calc_hp(Th, points); 
%
%  Parameters:  Th     - Pointer to the transmit aperture.  
%               points - Field points. Vector with three columns (x,y,z) 
%                        and one row for each field point. 
%
%  Return:      hp         - Emitted pressure field
%               start_time - The time for the first sample in field.
%
%  Version 1.0, November 28, 1995 by Joergen Arendt Jensen

function [hhp, start_time] = calc_hp (Th, points)

%  Check the point array

  [m,n]=size(points);
  if (n ~= 3)
    error ('Points array must have three columns');
    end

%  Call the C-part of the program to show aperture

  [hhp, start_time] = Mat_field (4002,Th,points);


