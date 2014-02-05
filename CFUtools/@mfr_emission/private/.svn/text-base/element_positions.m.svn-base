function obj = element_positions(obj)


%-----------------------------
% Element Position
%-----------------------------
if ~strcmp(obj.type,'custom')
    if obj.field_2_compat_mode
         % position axes (inclusive dead rows)
    obj.pos_x = (((-(obj.no_elm_x-1)/2):1:((obj.no_elm_x-1)/2))*obj.pitch_x)';
    obj.pos_y = (((-(obj.no_elm_y+obj.no_dead_rows-1)/2):((obj.no_elm_y+obj.no_dead_rows-1)/2))*obj.pitch_y)';
    % remove dead rows
    obj.pos_y(obj.dead_row_idx) = [];
    % count elm in the x-dir first
    idx_x = repmat((1:obj.no_elm_x)',obj.no_elm_y,1); %[1 2 3.. no_elm_x; 1 2.. no_elm_x ..]
    idx_y = reshape(ones(obj.no_elm_x,1)*(1:obj.no_elm_y),[],1); %[1 1 1.. 1; 2 2.. 2; 3 3.. 3; ...]
                                                                           % assemble position coords into [x y z] matrix          
    obj.pos = [obj.pos_x(idx_x), obj.pos_y(idx_y), zeros(obj.no_elm,1)];
    
    
    else  % use the new convention
        
        
    % position axes (inclusive dead rows)
    obj.pos_x = (((-(obj.no_elm_x+obj.no_dead_rows-1)/2):((obj.no_elm_x+obj.no_dead_rows-1)/2))*obj.pitch_x)';
    obj.pos_y = (((-(obj.no_elm_y-1)/2):1:((obj.no_elm_y-1)/2))*obj.pitch_y)';
    % remove dead rows
    obj.pos_x(obj.dead_row_idx) = [];
    % count elm in the x-dir first
    obj.pos_idx_x = repmat((1:obj.no_elm_x)',obj.no_elm_y,1); %[1 2 3.. no_elm_x; 1 2.. no_elm_x ..]
    obj.pos_idx_y = reshape(ones(obj.no_elm_x,1)*(1:obj.no_elm_y),[],1); %[1 1 1.. 1; 2 2.. 2; 3 3.. 3; ...]
    obj.pos = [obj.pos_x(obj.pos_idx_x), obj.pos_y(obj.pos_idx_y), zeros(obj.no_elm,1)];
    
    end
end

% if not set, set origin to center of transducer
if isempty(obj.origin_coord) && isempty(obj.origin_elm)
    obj.origin_coord = mean(obj.pos);
elseif ~isempty(obj.origin_elm)
    % origin
    dim1 = size(obj.origin_elm);
    if ~(max(dim1) == 2) && ~(min(dim1) == 1)
        error ('origin_elm must be a vector with two elm [idx_x idx_y].');
    end
    % from index to coord
    origin_x = interp1(obj.pos_x, obj.origin_elm(1));
    origin_y = interp1(obj.pos_y, obj.origin_elm(2));
    obj.origin_coord = [origin_x origin_y 0]; 
else
    dim1 = size(obj.origin_coord);
    if ~(max(dim1) == 3) && ~(min(dim1) == 1)
        error ('origin_coord must be a vector with three elm [x y z].');
    end                
end
% remove missing elm
% $$$ obj.pos(obj.elm_missing_idx,:) = [];
% $$$ obj.no_elm = size(obj.pos,1);

