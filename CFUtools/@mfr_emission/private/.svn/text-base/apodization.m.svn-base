function obj = apodization (obj)



trim_missing_elm_apo = false;





%% Circular Apodization (norm2)
if ~isempty(obj.no_act_elm)
    if obj.no_act_elm > obj.no_elm, 
        error('''no_act_elm'' cannot be larger than ''no_elm''');
    end 
    % determine distance from each element to apo origin
    dist = sqrt(sum((repmat(obj.origin_coord,[obj.no_elm 1])-obj.pos).^2,2));
    %Find number of active elm closest to apo_origin
    [dist_sort dist_sort_idx] = sort(dist);
    apo_idx = dist_sort_idx(1:obj.no_act_elm);
    chosen  = obj.pos(apo_idx, :);
    norm1   = dist_sort(1:obj.no_act_elm);
    % Get APO values
    index_norm =  norm1/max(norm1(:)); %[0-1]
    window     = mfr_window(obj.window, index_norm, obj.window_alpha);
    
    
    % Make sure minimum apo value is larger than threshhold
    if obj.window_thr > min(window)
        f=@(scale) window_min_val_residual(obj.window, index_norm, obj.window_alpha, ...
                                           obj.window_thr, scale);
        scale_0 = 1;
        [scale val exitflag output]  = fminsearch(f, scale_0);%, fmin_options)
        index_norm = index_norm*scale;
        window = mfr_window(obj.window, index_norm, obj.window_alpha);
    end   
    % Set max apo val = 1
    window = window/max(window(:));
    
    %output APO
    obj.apo = zeros(obj.no_elm,1);
    obj.apo(apo_idx(1:obj.no_act_elm)) = window(1:obj.no_act_elm);

    trim_missing_elm_apo = false;

    
else % obj.no_act_elm is empty
    % Make sure number of active elm is sane
    if (isempty(obj.no_act_elm_x) || isempty(obj.no_act_elm_y))
        error(['''no_act_elm'' or ''no_act_elm_x'' and ''no_act_elm_y'' must be set when not using the ' ...
               '''custom'' array.']); 
    end
    if obj.no_act_elm_y > obj.no_elm_y, 
        error('''no_act_elm_y'' is larger than ''no_elm_y''');
    end 
    if obj.no_act_elm_x > obj.no_elm_x, 
        error('''no_act_elm_x'' is larger than ''no_elm_x''');
    end 
    
    
    
    %% normal transducers
    if ~isempty(obj.apo_x) && ~isempty(obj.apo_y)
        if length(obj.apo_x) ~= obj.no_elm_x || length(obj.apo_x) ~= obj.no_elm_x
            error(['''apo_x'' and ''apo_y'' must have the same length as ''no_elm_x''  ' ...
                   'respectively ''no_elm_y'''])
        end
        obj.apo = obj.apo_x *  obj.apo_y';
        obj.apo = obj.apo(:);
        trim_missing_elm_apo = true;

    elseif isempty(obj.apo)
        pos_x_tmp = obj.pos_x;
        pos_y_tmp = obj.pos_y;
        
        % find nearest elm to origin -- X-dim
        obj.apo_x_idx = zeros(1,obj.no_act_elm_x);
        %TODO: rewrite using sort!
        for idx = 1:obj.no_act_elm_x
            % find nearest element in tmp array
            [min_dist min_idx] = min (abs(pos_x_tmp - obj.origin_coord(1)));  %#ok
                                                                              % find corresponding element index
            elm_idx = find (obj.pos_x == pos_x_tmp(min_idx));
            % save element index
            obj.apo_x_idx(idx) = elm_idx;
            % remove found element from future searches
            pos_x_tmp(min_idx) = [];
        end
        
        % find neares elm to origin -- Y-dim
        obj.apo_y_idx = zeros(1,obj.no_act_elm_y);
        %TODO: rewrite using sort!
        for idx = 1:obj.no_act_elm_y
            % find nearest element in tmp array
            [min_dist min_idx] = min (abs(pos_y_tmp - obj.origin_coord(2))); %#ok
                                                                             % find corresponding element index
            elm_idx = find (obj.pos_y == pos_y_tmp(min_idx));
            % save element index
            obj.apo_y_idx(idx) = elm_idx;
            % remove found element from future searches
            pos_y_tmp(min_idx) = [];
        end
        % Sort array
        obj.apo_x_idx = sort(obj.apo_x_idx);
        obj.apo_y_idx = sort(obj.apo_y_idx);
        
        % Index norm [0-1]
        index_norm_x = obj.pos_x(obj.apo_x_idx) - min(obj.pos_x(obj.apo_x_idx));
        if max(index_norm_x) > 0
            index_norm_x = abs(2*index_norm_x/max(index_norm_x)-1); % inverse triangle function ex: [1, 0.5, 0, 0.5, 1]
        end
        index_norm_y = obj.pos_y(obj.apo_y_idx) - min(obj.pos_y(obj.apo_y_idx));
        if max(index_norm_y) > 0
            index_norm_y = abs(2*index_norm_y/max(index_norm_y)-1); % inverse triangle function ex: [1, 0.5, 0, 0.5, 1]
        end
        
        
        
        % Window values. Make sure minimum APO value is larger or equal to threshold
        window_x = mfr_window(obj.window, index_norm_x, obj.window_alpha);
        if obj.window_thr > min(window_x)
            f=@(scale) window_min_val_residual(obj.window, index_norm_x, obj.window_alpha, obj.window_thr,scale);
            scale_0 = 1;
            scale   = fminsearch(f, scale_0);
            index_norm_x = index_norm_x*scale;
            window_x = mfr_window(obj.window, index_norm_x, obj.window_alpha);
        end

        window_y = mfr_window(obj.window, index_norm_y, obj.window_alpha);
        if obj.window_thr > min(window_y)
            f=@(scale) window_min_val_residual(obj.window, index_norm_y, obj.window_alpha, obj.window_thr,scale);
            scale_0 = 1;
            scale   = fminsearch(f, scale_0);
            index_norm_y = index_norm_y*scale;
            window_y = mfr_window(obj.window, index_norm_y, obj.window_alpha);
        end

        % Fill APO with zeros
        obj.apo_x = zeros(1,obj.no_elm_x);
        obj.apo_y = zeros(1,obj.no_elm_y);
        % Write active window to APO
        obj.apo_x(obj.apo_x_idx) = window_x;
        obj.apo_y(obj.apo_y_idx) = window_y;
        % Make sure maximum is 1
        obj.apo_x = obj.apo_x/max(obj.apo_x(:));
        obj.apo_y = obj.apo_y/max(obj.apo_y(:));
        % vectors to matrix
        obj.apo = obj.apo_x'*obj.apo_y;
        % matrix to vector
        obj.apo = obj.apo(:)/max(obj.apo(:));

        
    else % obj.apo is not empty
        if length(obj.apo) == (size(obj.pos,1))
            trim_missing_elm_apo = false; %input apo already has the right size.
        elseif length(obj.apo) == (size(obj.pos,1) + length(obj.missing_elm_idx))
            trim_missing_elm_apo = true; %input apo must be trimmed
        else 
            error ('Apodization must have same length as number of elm.')
        end
    end
    
    % make sure max of apo is one
    if max(obj.apo) ~= 0
        obj.apo = obj.apo/max(obj.apo);
    end
end









% remove any missing/dead/inactive elements
if ~isempty([obj.elm_missing_idx obj.elm_dead_idx]) && trim_missing_elm_apo
    obj.apo([obj.elm_missing_idx obj.elm_dead_idx])   = [];
end




