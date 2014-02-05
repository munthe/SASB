function obj = aperture_setup(obj)

if (strcmp(obj.type,'vermon_32x32')) % Full SARUS
    obj.no_elm_x = 32;
    obj.no_elm_y = 32;
    obj.no_elm   = obj.no_elm_x * obj.no_elm_y; 
    obj.no_dead_rows = 3;
    obj.dead_row_idx = [27 18 9];
    
    
    % ---------------------------------------------
    % Simulation Probes
    % ---------------------------------------------
elseif (strcmp(obj.type,'dense_custom'))
    if isempty(obj.no_elm_x) || isempty(obj.no_elm_y)
        error('''no_elm_x'' and ''no_elm_y'' must be set when using the ''dense_custom'' array.');
    end
    obj.no_elm = obj.no_elm_x * obj.no_elm_y; 
    obj.no_dead_rows = 0;
    obj.dead_row_idx = [];
    
elseif (strcmp(obj.type,'dense_32x32'))
    obj.no_elm_x = 32;
    obj.no_elm_y = 32;
    obj.no_elm = obj.no_elm_x * obj.no_elm_y; 
    obj.no_dead_rows = 0;
    obj.dead_row_idx = [];
    
elseif (strcmp(obj.type,'dense_64x64'))
    obj.no_elm_x = 64;
    obj.no_elm_y = 64;
    obj.no_elm = obj.no_elm_x * obj.no_elm_y; 
    obj.no_dead_rows = 0;
    obj.dead_row_idx = [];

    
    
    
    % ---------------------------------------------
    % Custom Probe
    % ---------------------------------------------
elseif (strcmp(obj.type,'custom'))
    %Everything determined externally
    if isempty(obj.pos)
        error('By ''custom'' transducer the element position (pos) must be given.')
    end
    if isempty(obj.no_act_elm) && isempty(obj.apo)
        error('When using ''custom'' transducer and no predefined apodization, the number of active elm (no_act_elm) must be given.')
    end
    obj.no_elm = size(obj.pos,1);
    
else
    error([' Unsupported aperture type.\n Known apertures are: "custom", "vermon_32x32", ' ...
           '"dense_32x32"," "dense_64x64 and "dense_custom".']);
end


