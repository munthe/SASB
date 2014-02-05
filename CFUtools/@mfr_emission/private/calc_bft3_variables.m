function obj = calc_bft3_variables(obj)
% This function calculates values that eases the task of beamforming with BFT3.

% Index of element with delay==0 and apodization > 0
cf_idx = find((obj.delays==0).*(obj.apo), 1,'first');


%% VS behind aperture
if (obj.focus_r <= 0)
    obj.bft3_focus        = [];
    obj.bft3_center_focus = obj.vs;
    obj.bft3_tstart_delta = norm(obj.pos(cf_idx,:) - obj.vs)/obj.c;
    return
end




%% VS in front of aperture
% the naive
if obj.bft3_time_comp_mode == 1  
    obj.bft3_focus        = obj.vs;
    obj.bft3_center_focus = obj.pos(cf_idx,:);
    obj.bft3_tstart_delta = 0;

    % The advanced
elseif obj.bft3_time_comp_mode ==2
    obj.bft3_focus        = obj.vs;
    obj.bft3_center_focus = obj.origin_coord;
    obj.bft3_tstart_delta = abs(norm( obj.origin_coord  -obj.vs ) - ...
                                norm( obj.pos(cf_idx,:) -obj.vs ))/obj.c;
    
    % When VS = CF (single element emissions)
elseif obj.bft3_time_comp_mode ==3
    obj.bft3_focus        = [];
    obj.bft3_center_focus = obj.vs;
    obj.bft3_tstart_delta = norm(obj.pos(cf_idx,:) - obj.vs)/obj.c;

else 
    error('bft3_time_comp_mode must be 1,2 or 3.')
end
