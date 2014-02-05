function residual = window_min_val_residual(window, index_norm, window_alpha, goal_val,scale)

window   = mfr_window(window, index_norm*scale, window_alpha);
min_val  = min(window);
max_val = max(window);
residual = (min_val - goal_val)^2;

%fprintf('scale= %f\nmin  = %f\nmax  = %f\n\n', scale, min_val, max_val)


