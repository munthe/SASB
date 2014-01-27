%  Compress the data to show 60 dB of
%  dynamic range
%
%  Example by Joergen Arendt Jensen, Nov. 28, 1995.

%  Adjust the data in time and display it as
%  a gray scale image
%
%  Version 1.1, April 1, 1998, Joergen Arendt Jensen
%  Version 1.2, August 13, 2007, JAJ: 
%       Printout changed
%       Display to a 60 dB dynamic range

disp('Making display of image')
D=40;             %  Sampling frequency decimation factor
min_sample=9/1000/c;
max_sample=max([max(times)*fs 2*0.121/c*fs]);
[n,m]=size(image_data);
n=n+(max_sample-min_sample);
env=zeros(round(n/D),no_lines);
for i=1:no_lines
  rf_env=abs(hilbert([zeros(round(times(i)*fs-min_sample),1); image_data(:,i)]));
  rf_env=rf_env(1:D:max(size(rf_env)));
  env(1:max(size(rf_env)),i)=rf_env;
  end

%  Do logarithmic compression to a 60 dB dynamic range

log_env=20*log10(env);
log_env=log_env-max(max(log_env)) + 60;
log_env=127*log_env/60;

%  Set values -80 dB to -80 dB

log_env=(log_env > -80).*log_env + (log_env < -80)*(-80);

%  Make an interpolated image

ID=5;
[n,m]=size(log_env);
new_env=zeros(n,m*ID);
for i=1:n
  new_env(i,:)=interp(log_env(i,:),ID);
  end
[n,m]=size(new_env);
  
fn=fs/D;
image(((1:(ID*no_lines-1))*d_x/ID-no_lines*d_x/2)*1000,((1:n)/fn+min_sample/fs)*1540/2*1000,new_env)
colormap(gray(128))
axis('image')
axis([-10 10 10 120])
axis off
drawnow
