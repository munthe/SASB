%  Function for making the interpolation of an ultrasound image.
%  The routine make_tables must have been called previously.
%
%  Input parameters: 
%
%         envelope_data - The envelope detected and log-compressed data as
%                         an integer array as 8 bits values  (uint8)
%
%  Output:  img_data    - The final image as 8 bits values
%
%  Calling: img_data = make_interpolation (envelope_data); 
%
%  Version 1.0, 14/2-1999, JAJ
%  Version 1.1,  8/9-2003, JAJ: Error in help text corrected.
%  Version 1.2, 29/8-2011, JAJ: Unsigned 32 bits version

function img_data = make_interpolation (envelope_data)

%  Call the appropriate function


img_data = fast_int_uint32 (2, envelope_data);
	    
