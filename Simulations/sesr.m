% Conventional imaging with single emision single receive focus and
% apodixation vector
%

x= -image_width/2;
image_data=0;

for i=1:no_lines
   % Set the focus for this direction

  xdc_center_focus (emit_aperture,[x 0 0]);
  xdc_focus (emit_aperture, 0, [x 0 trans_focus]);
  xdc_center_focus (receive_aperture,[x 0 0]);
  xdc_focus (receive_aperture, 0, [x 0 receive_focus]);

   % Calculate the apodization 
   
  N_pre  = round(x/(width+kerf) + N_elements/2 - N_active/2);
  N_post = N_elements - N_pre - N_active;
  apo_vector=[zeros(1,N_pre) apo zeros(1,N_post)];
  xdc_apodization (emit_aperture, 0, apo_vector);
  xdc_apodization (receive_aperture, 0, apo_vector);

  % Calculate the received response
  [v, t1]=calc_scat(emit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);

  % Store the result
  image_data(1:max(size(v)),i)=v;
  times(i) = t1;

  % Steer in another angle
  x = x + d_x;
  
end
