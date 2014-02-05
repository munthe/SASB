function obj = remove_dead_elements(obj)

% Get missing elements from missing DAUP index
if ~isempty(obj.daup_missing_idx)
    elm_missing_idx = daup_idx_to_ch_idx(obj.daup_missing_idx);
    obj.elm_missing_idx = [obj.elm_missing_idx(:); elm_missing_idx(:)];
end

% make sure only one of each index
obj.elm_missing_idx = unique(obj.elm_missing_idx);

% remove position and delays
obj.pos(obj.elm_missing_idx,:)    = [];
obj.no_elm                        = size(obj.pos,1);
%obj.delays(obj.elm_missing_idx,:) = [];

% remove apo
if obj.apo_external_set
    obj.apo(obj.elm_missing_idx) = [];
end
