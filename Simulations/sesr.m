% Conventional imaging with single emision single receive focus and
% apodixation vector
%
% Calling: 
%
% Output: 
%
% 
% function [img] = sesr(phantom,emis_focus,trans_focus,apo_vector)
% 
% 
% 
% end


%  Make imaging of phantom data for a single transmit/receive focus
%
%  Version 1.2, August 14, 1998, Joergen Arendt Jensen
%       Problem with focusing reference fixed

%  Make the apodization vector

apo=ones(1,N_active);

z_focus=60/1000;          %  Transmit focus

x= -image_width/2;
image_data=0;
for i=1:no_lines
   %   Set the focus for this direction

  xdc_center_focus (emit_aperture,[x 0 0]);
  xdc_focus (emit_aperture, 0, [x 0 z_focus]);
  xdc_center_focus (receive_aperture,[x 0 0]);
  xdc_focus (receive_aperture, 0, [x 0 z_focus]);

   %  Calculate the apodization 
   
  N_pre  = round(x/(width+kerf) + N_elements/2 - N_active/2);
  N_post = N_elements - N_pre - N_active;
  apo_vector=[zeros(1,N_pre) apo zeros(1,N_post)];
  xdc_apodization (emit_aperture, 0, apo_vector);
  xdc_apodization (receive_aperture, 0, apo_vector);

  %   Calculate the received response

  [v, t1]=calc_scat(emit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);

  %  Store the result

  image_data(1:max(size(v)),i)=v;
  times(i) = t1;

  %  Steer in another angle

  x = x + d_x;
  end
