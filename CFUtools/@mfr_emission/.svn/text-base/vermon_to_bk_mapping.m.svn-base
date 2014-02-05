% This file permutates the input vector. It changes the element
% ordering from the Vermon standard to BK-Med. standard.
% Vermon standard: Transducer mark, marks the direction of positive y
% BK-Med standard: Transducer mark, marks the direction of positive x

function ch_idx = elm_idx_to_ch_idx(elm_idx)



ch_map = 
obj.pos_x = (((-(obj.no_elements_x+obj.no_dead_rows-1)/2):((obj.no_elements_x+obj.no_dead_rows-1)/2))*obj.pitch_x)';
obj.pos_y = (((-(obj.no_elements_y-1)/2):1:((obj.no_elements_y-1)/2))*obj.pitch_y)';
% remove dead rows
obj.pos_x(obj.dead_row_idx) = [];
% count elements in the x-dir first
idx_x = repmat((1:obj.no_elements_x)',obj.no_elements_y,1); %[1 2 3.. no_elements_x; 1 2.. no_elements_x ..]
idx_y = reshape(ones(obj.no_elements_x,1)*(1:obj.no_elements_y),[],1); %[1 1 1.. 1; 2 2.. 2; 3 3.. 3; ...]


ch_idx = 0;









