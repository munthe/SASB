%  Function for calculating the different weight tables
%  for the interpolation. The tables are calculated according
%  to the parameters in this function and are stored in
%  the internal C-code.
%
%  Input parameters: 
%
%        start_depth    - Depth for start of image in meters
%        image_size     - Size of image in meters
%		   
%	 start_of_data  - Depth for start of data in meters
%	 delta_r        - Sampling interval for data in meters
%	 N_samples      - Number of samples in one envelope line
%		   
%	 theta_start    - Angle for first line in image
%	 delta_theta    - Angle between individual lines
%	 N_lines        - Number of acquired lines
%		   
%	 scaling        - Scaling factor from envelope to image
%	 Nz             - Size of image in pixels
%	 Nx             - Size of image in pixels
%
%  Output:  Nothing, everything is stored in the C program
%
%  Calling: make_tables (start_depth, image_size,             ...
%		         start_of_data, delta_r, N_samples,   ...
%		         theta_start, delta_theta, N_lines,   ...
%		         scaling, Nz, Nx); 
%
%  Version 1.0, 14/2-1999, JAJ
%  Version 1.1, 11/9-2001, JAJ: Help text corrected

function res = make_tables (start_depth, image_size,             ...
                            start_of_data, delta_r, N_samples,   ...
		            theta_start, delta_theta, N_lines,   ...
		            scaling, Nz, Nx)

%  Call the appropriate function


fast_int (1,start_depth, image_size,             ...
            start_of_data, delta_r, N_samples,   ...
	    theta_start, delta_theta, N_lines,   ...
	    scaling, Nz, Nx);
	    
%  Return nothing

ret=0;

