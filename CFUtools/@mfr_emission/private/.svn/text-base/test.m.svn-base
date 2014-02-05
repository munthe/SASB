function obj = test(obj)

methods(obj);
    properties(obj);

test_obj = mfr_emission('c',1540, ...
                        'no_act_elm',1024,...
                        'apo',[zeros(200,1);ones(195,1);zeros(603,1)],...
                        'focus_r', 60/1000, ...
                        'zx_axis', 30*pi/180, ...
                        'zy_axis', 30*pi/180,...
                        'elm_missing_idx',[210:214 (700:720)]);

figure
test_obj.plot_vs;

figure;
plot(test_obj.delays)

help mfr_emission

    
    
    
    