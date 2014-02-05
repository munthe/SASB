%% This function calculates the delay of each transducer element.

function obj = calc_delays(obj)


% $$$ % force angle to be [0 2*pi[  (0 is straight down)
% $$$ zy_eff = mod(obj.zy_axis,2*pi);
% $$$ zx_eff = mod(obj.zx_axis,2*pi);            
% $$$ 
% $$$ % force angle to be ]-pi pi]
% $$$ if zy_eff > pi, zy_eff = 2*pi - zy_eff; end;
% $$$ if zx_eff > pi, zx_eff = 2*pi - zx_eff; end;

if abs(obj.zy_axis) > pi/2  ||  abs(obj.zx_axis) > pi/2
    % Too large angles
    error(['   Virtual source behind transducer implemented as negativ focus_r - not as steering ' ...
           'angles behind aperture  (zy=%.2f zx=%.2f)']', ...
          zy_eff*180/pi, zx_eff*180/pi);
end    



% Calc time-of-flight from element to focus point
obj.delays = sqrt(sum((obj.pos - repmat(obj.vs,[obj.no_elm 1])).^2,2))/obj.c;

% If virtual source in front of the transducer
if obj.focus_r >= 0; 
    obj.delays = -obj.delays;    
end
% If virtual source behind the transducer -> do nothing, delays are
% already correct.


% remove delays from non-transmitting elements
obj.delays(obj.apo==0) = min(obj.delays(obj.apo>0));
obj.delays =  obj.delays - min(obj.delays);

