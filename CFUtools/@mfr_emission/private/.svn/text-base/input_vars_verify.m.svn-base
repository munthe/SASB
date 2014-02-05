% This script verifys a sane state of the input variables


% set speed of sound, if water_temp is set.
if (isempty(obj.c) && ~isempty(obj.water_temp))
    obj.c = sound_speed_calc(obj.water_temp);
end

% Sound speed test
if isempty(obj.c)
    error([' You must set the speed of sound.\n If measuring phantoms in water, ', ...
           'it is recommended to set the variable "water_temp". This will then automatically ', ...
           'set the correct speed of sound in water.'])
end


if obj.window_thr < 0 || obj.window_thr >= 1, 
    error('window_thr must be larger than or equal to 0, and smaller than 1.');
end
if ~isempty(st.no_act_elm) && (~isempty(st.no_act_elm_x) || ~isempty(st.no_act_elm_y))
    error('When ''no_act_elm'' is set, ''no_act_elm_x'' and ''no_act_elm_y'' cannot be set.');
end




if (strcmp(obj.type,'dense_custom'))
    if isempty(obj.no_elements_x) || isempty(obj.no_elements_x)
        error(['''no_elements_x'' and ''obj.no_elements_y'' must be set when using the ''dense_custom'' ' ...
               'aperture type']);
    end
end

