


function [angles] = cfu_setup_3d_lines (varargin)

%
% Sets up the angle of each transmit and receive line. (Mostly for explososcan)
% 
% [angles angles_tx angles_rx] = cfu_setup_3d_lines (varargin)
%
% 
% -- Input Arguments --
% 'span_zy': span (interval) of angles of rotation around x-axis [rad]
% 'span_zx': span (interval) of angles of rotation around y-axis [rad]
% 'no_TX_zy': Number og TX lines in zy dimension.
% 'no_TX_zx': Number og TX lines in zx dimension.
% 'no_lines_per_tx_zy': Number of RX lines per TX line.
% 'no_lines_per_tx_zx': Number of RX lines per TX line.
% 'no_overlap':  Number of lines that should overlap from two different
% -- Output --
% Per Transmission Info:
% 'angles{idx}.tx_zy': X-rotation for emission idx.
% 'angles{idx}.tx_zx': Y-rotation for emission idx.
% 'angles{idx}.rx_zy': Array of X-rotations for receive lines centered on tx_zy.
% 'angles{idx}.rx_zx': Array of Y-rotations for receive lines centered on tx_zx.
% 
% Entire image info:
% 'angles{idx}.rx_zy_ar': Array of all receive  rot_zy angles for entire image
% 'angles{idx}.rx_zx_ar': Array of all receive  rot_zx angles for entire image
% 'angles{idx}.tx_zy_ar': Array of all transmit rot_zy angles for entire image
% 'angles{idx}.tx_zx_ar': Array of all transmit rot_zx angles for entire image
%
% MFR 2011-10-02 - Init version.
% MFR 2011-10-18 - Entire image info added.
% MFR 2011-11-28 - Overlapping lines and weight of these added. rot_x renamed to zy.
% MFR 2012-06-05 - Now also works in 1D (without overlap).
% MFR 2013-08-05 - Renamed from mfr_explosolines_setup to cfu_setup_3d_lines.
%

%
% OBS: keep 

%
% A note of notation: in this file the following notation is often used:
% rx_x which is short for rx_rot_x -- in other words rotation around the x-axis in radians.
% This means in cartesian coordinates the y-coordinate changes when rx_x is changed. 
%


st.span_zy            = [];
st.span_zx            = [];
st.no_TX_zy           = 1;
st.no_TX_zx           = 1;
st.no_lines_per_tx_zy = 1;
st.no_lines_per_tx_zx = 1;
st.no_overlap         = 0;
% get input
st = cfu_parse_input_parameters(st,varargin);
%st = bft3_va_arg(st,varargin);

% remove struct
span_zy            = st.span_zy;
span_zx            = st.span_zx;
no_TX_zy           = st.no_TX_zy;
no_TX_zx           = st.no_TX_zx;
no_lines_per_tx_zy = st.no_lines_per_tx_zy;
no_lines_per_tx_zx = st.no_lines_per_tx_zx;
no_overlap         = st.no_overlap;

%% Rotation X RX array
if length(span_zy) ==1
    no_lines_zy = 1;
    rx_rot_zy   = span_zy;
    if (no_TX_zy > 1), 
        error(sprintf(['no_TX_zy >1, but no span of angles is given.\n', ...
                      'Either set no_TX_zy=1 or define an angle span (span_zy).'])); 
    end
else
    % interval info
    no_lines_zy = no_TX_zy*no_lines_per_tx_zy;
    int_zy_width  = abs(span_zy(2)-span_zy(1));
    int_zy_center = mean(span_zy);
    int_zy_rising = sign(span_zy(2)-span_zy(1)); %normal increasing interval
    % fill receive array
    if int_zy_rising
        if no_lines_zy ~=1
            rx_rot_zy = (1:no_lines_zy)/(no_lines_zy-1);
        else
            rx_rot_zy = (1:no_lines_zy);
        end
    else
        if no_lines_zy ~=1
            rx_rot_zy = (no_lines_zy:-1:1)/(no_lines_zy-1);
        else
            rx_rot_zy = (no_lines_zy:-1:1);
        end
    end
    % normalise width to "_width" and center array at "_center"
    rx_rot_zy = (rx_rot_zy - mean(rx_rot_zy))*int_zy_width + int_zy_center;
end


%% Rotation Y RX array
if length(span_zx) ==1
    no_lines_zx = 1;
    rx_rot_zx   = span_zx;
    if (no_TX_zx > 1), 
        error(sprintf(['  no_TX_zx >1, but no span of angles is given.\n', ...
                      '  Either set no_TX_zx=1 or define an angle span (span_zx).'])); 
    end
else
    no_lines_zx   = no_TX_zx*no_lines_per_tx_zx;
    int_zx_width  = abs(span_zx(2)-span_zx(1));
    int_zx_center = mean(span_zx);
    int_zx_rising = sign(span_zx(2)-span_zx(1)); %normal increasing interval

    % fill receive array
    if int_zx_rising
        rx_rot_zx = (1:no_lines_zx)/(no_lines_zx-1);
    else
        rx_rot_zx = (no_lines_zx:-1:1)/(no_lines_zx-1);
    end
    
    % normalise width to "_width" and center array at "_center"
    rx_rot_zx = (rx_rot_zx - mean(rx_rot_zx))*int_zx_width + int_zx_center;
end




%% X-rot TX array
% Transmit array sample points 
start_zy = (no_lines_per_tx_zy-1)/2 +1; %ceil(no_lines_per_tx_zy/2);
dy       = no_lines_per_tx_zy;
end_zy   = start_zy + dy*(no_TX_zy-1);
tx_zy_sample_p = (start_zy:dy:end_zy); %sample point
if no_lines_zy == 1 %special case (interp1 needs atleast 2 points)
    tx_rot_zy = rx_rot_zy(tx_zy_sample_p);
else
    % Interpolate RX array
    tx_rot_zy = interp1(rx_rot_zy,tx_zy_sample_p);
end


%% Y-rot TX array
% Transmit array sample points 
start_zx = (no_lines_per_tx_zx-1)/2 +1;
dx      = no_lines_per_tx_zx;
end_zx   = start_zx + dx*(no_TX_zx-1);
tx_zx_sample_p = (start_zx:dx:end_zx);
if no_lines_zx == 1 %special case (interp1 needs atleast 2 points)
    tx_rot_zx = rx_rot_zx(tx_zx_sample_p);
else
    % Interpolate RX array
    tx_rot_zx = interp1(rx_rot_zx,tx_zx_sample_p);
end






% Cluster rx lines into groups centered around each TX line
% X
for tx_idx = 1:length(tx_rot_zy)
    rx_zy_idx{tx_idx}  = (1:no_lines_per_tx_zy) + (tx_idx-1)*no_lines_per_tx_zy;
    % add overlapping lines (index)
    rx_zy_idx{tx_idx}  = [rx_zy_idx{tx_idx}(1)-(no_overlap:-1:1) rx_zy_idx{tx_idx} rx_zy_idx{tx_idx}(end)+(1:no_overlap)];
    rx_zy_overlap{tx_idx} = [ones(1,no_overlap) zeros(1,no_lines_per_tx_zy) ones(1,no_overlap)];
    % remove out of bounds overlapping lines
    %fprintf('x_idx pre %i ',length(rx_zy_idx{tx_idx}));
    remove_idx = rx_zy_idx{tx_idx}<1 | rx_zy_idx{tx_idx}>length(rx_rot_zy);
    rx_zy_idx{tx_idx}(remove_idx) = [];
    rx_zy_overlap{tx_idx}(remove_idx) = [];
    %fprintf(' post: %i\n',length(rx_zy_overlap{tx_idx}));
    % pick the angles (line idx to line angle)
    rx_zy{tx_idx}      = rx_rot_zy(rx_zy_idx{tx_idx});
    % distance from TX to RX in angles
    rx_zy_dist{tx_idx} = abs(tx_rot_zy(tx_idx) - rx_zy{tx_idx});
end

% Y
for tx_idx = 1:length(tx_rot_zx)
    rx_zx_idx{tx_idx} = (1:no_lines_per_tx_zx) + (tx_idx-1)*no_lines_per_tx_zx;
    % add overlapping lines (index)
    rx_zx_idx{tx_idx}     = [rx_zx_idx{tx_idx}(1)-(no_overlap:-1:1) rx_zx_idx{tx_idx} rx_zx_idx{tx_idx}(end)+(1:no_overlap)];
    rx_zx_overlap{tx_idx} = [ones(1,no_overlap) zeros(1,no_lines_per_tx_zx) ones(1,no_overlap)];
    % remove out of bounds overlapping lines
    remove_idx = rx_zx_idx{tx_idx}<1 | rx_zx_idx{tx_idx}>length(rx_rot_zx);
    rx_zx_idx{tx_idx}(remove_idx)     = [];
    rx_zx_overlap{tx_idx}(remove_idx) = [];
    % pick the angles (line idx to line angle)
    rx_zx{tx_idx}     = rx_rot_zx(rx_zx_idx{tx_idx});
    % distance from TX to RX in angles
    rx_zx_dist{tx_idx} = abs(tx_rot_zx(tx_idx) - rx_zx{tx_idx});
end




if st.no_overlap ~0
    %% The weight of each line
    % Since the lines are distributed equaly, the distances are not calculated in angles, but in indices.

    %Distance between sources (in indices)
    dist_to_next_source_zy = no_lines_per_tx_zy;
    dist_to_next_source_zx = no_lines_per_tx_zx;
    %maximal distance from source to line that should be weighted
    %max_dist_zy = floor(no_lines_per_tx_zy/2) + no_overlap;
    %max_dist_zx = floor(no_lines_per_tx_zx/2) + no_overlap;

    % distance from source to line, within each group (lines to be beamformed per TX)
    % center column
    dist_src_5_zy = rx_zy_idx{2} - tx_zy_sample_p(2); %center source
    dist_src_4_zy = dist_src_5_zy;
    dist_src_6_zy = dist_src_5_zy;
    % right column
    dist_src_7_zy = dist_to_next_source_zy - dist_src_5_zy;
    dist_src_8_zy = dist_to_next_source_zy - dist_src_5_zy;
    dist_src_9_zy = dist_to_next_source_zy - dist_src_5_zy;
    % left column
    dist_src_1_zy = dist_to_next_source_zy + dist_src_5_zy;
    dist_src_2_zy = dist_to_next_source_zy + dist_src_5_zy;
    dist_src_3_zy = dist_to_next_source_zy + dist_src_5_zy;

    % center row
    dist_src_5_zx = rx_zx_idx{2} - tx_zx_sample_p(2); %center source
    dist_src_2_zx = dist_src_5_zx;
    dist_src_8_zx = dist_src_5_zx;
    % bottom row
    dist_src_3_zx = dist_to_next_source_zx - dist_src_5_zx;
    dist_src_6_zx = dist_to_next_source_zx - dist_src_5_zx;
    dist_src_9_zx = dist_to_next_source_zx - dist_src_5_zx;
    % top column
    dist_src_1_zx = dist_to_next_source_zx + dist_src_5_zx;
    dist_src_4_zx = dist_to_next_source_zx + dist_src_5_zx;
    dist_src_7_zx = dist_to_next_source_zx + dist_src_5_zx;

    % vectors to norm 2 matrices
    dist_src_mtr(:,:,1) = sqrt(repmat((dist_src_1_zy.^2),   length(dist_src_1_zx),1) + ...
                               repmat((dist_src_1_zx.^2)',1,length(dist_src_1_zy)));
    dist_src_mtr(:,:,2) = sqrt(repmat((dist_src_2_zy.^2),   length(dist_src_2_zx),1) + ...
                               repmat((dist_src_2_zx.^2)',1,length(dist_src_2_zy)));
    dist_src_mtr(:,:,3) = sqrt(repmat((dist_src_3_zy.^2),   length(dist_src_3_zx),1) + ...
                               repmat((dist_src_3_zx.^2)',1,length(dist_src_3_zy)));
    dist_src_mtr(:,:,4) = sqrt(repmat((dist_src_4_zy.^2),   length(dist_src_4_zx),1) + ...
                               repmat((dist_src_4_zx.^2)',1,length(dist_src_4_zy)));
    dist_src_mtr(:,:,5) = sqrt(repmat((dist_src_5_zy.^2),   length(dist_src_5_zx),1) + ...
                               repmat((dist_src_5_zx.^2)',1,length(dist_src_5_zy)));
    dist_src_mtr(:,:,6) = sqrt(repmat((dist_src_6_zy.^2),   length(dist_src_6_zx),1) + ...
                               repmat((dist_src_6_zx.^2)',1,length(dist_src_6_zy)));
    dist_src_mtr(:,:,7) = sqrt(repmat((dist_src_7_zy.^2),   length(dist_src_7_zx),1) + ...
                               repmat((dist_src_7_zx.^2)',1,length(dist_src_7_zy)));
    dist_src_mtr(:,:,8) = sqrt(repmat((dist_src_8_zy.^2),   length(dist_src_8_zx),1) + ...
                               repmat((dist_src_8_zx.^2)',1,length(dist_src_8_zy)));
    dist_src_mtr(:,:,9) = sqrt(repmat((dist_src_9_zy.^2),   length(dist_src_9_zx),1) + ...
                               repmat((dist_src_9_zx.^2)',1,length(dist_src_9_zy)));



% $$$ max_val = max(max(dist_src_mtr(:,:,1)));
% $$$ figure;
% $$$ subplot(3,3,1);imagesc(dist_src_mtr(:,:,1), [0 max_val])
% $$$ subplot(3,3,4);imagesc(dist_src_mtr(:,:,2), [0 max_val])
% $$$ subplot(3,3,7);imagesc(dist_src_mtr(:,:,3), [0 max_val])
% $$$ subplot(3,3,2);imagesc(dist_src_mtr(:,:,4), [0 max_val])
% $$$ subplot(3,3,5);imagesc(dist_src_mtr(:,:,5), [0 max_val])
% $$$ subplot(3,3,8);imagesc(dist_src_mtr(:,:,6), [0 max_val])
% $$$ subplot(3,3,3);imagesc(dist_src_mtr(:,:,7), [0 max_val])
% $$$ subplot(3,3,6);imagesc(dist_src_mtr(:,:,8), [0 max_val])
% $$$ subplot(3,3,9);imagesc(dist_src_mtr(:,:,9), [0 max_val])



    %% number of sources contributing to each line within the group

    % spacial cases: corners, edges, and middle
    for idx_tx_zy = [1 2 no_TX_zy]
        for idx_tx_zx = [1 2 no_TX_zx]
            weight_tmp = working_loop (idx_tx_zy, idx_tx_zx, ...
                                       no_TX_zx, no_TX_zy, ...
                                       no_lines_per_tx_zx, no_lines_per_tx_zy,...
                                       no_overlap, dist_src_mtr);
            weight{idx_tx_zy,idx_tx_zx} = weight_tmp;
        end
    end

    % top edge
    idx_tx_zy = 1;
    for idx_tx_zx = 3:no_TX_zx-1
        weight{idx_tx_zy,idx_tx_zx} =  weight{1,2};
    end
    % bottom edge
    idx_tx_zy = no_TX_zy;
    for idx_tx_zx = 3:no_TX_zx-1
        weight{idx_tx_zy,idx_tx_zx} =  weight{no_TX_zy,2};
    end
    % left edge
    idx_tx_zx = 1;
    for idx_tx_zy = 3:no_TX_zy-1
        weight{idx_tx_zy,idx_tx_zx} =  weight{2,1};
    end
    % right edge
    idx_tx_zx = no_TX_zx;
    for idx_tx_zy = 3:no_TX_zy-1
        weight{idx_tx_zy,idx_tx_zx} =  weight{2,no_TX_zx};
    end
    % middle
    for idx_tx_zy = 2:no_TX_zy-1
        for idx_tx_zx = 2:no_TX_zx-1
            weight{idx_tx_zy,idx_tx_zx} = weight{2,2};
        end
    end

% $$$ figure;imagesc(weight{2,2})


    %% Create output cell array
    idx=0;
    for idx_zy = 1:no_TX_zy
        for idx_zx = 1:no_TX_zx
            idx = idx+1;
            angles{idx}.rx_zx     = rx_zx{idx_zx};
            angles{idx}.rx_zy     = rx_zy{idx_zy};
            angles{idx}.rx_zx_idx = rx_zx_idx{idx_zx};
            angles{idx}.rx_zy_idx = rx_zy_idx{idx_zy};
            angles{idx}.tx_zy     = tx_rot_zy(idx_zy);
            angles{idx}.tx_zx     = tx_rot_zx(idx_zx);
            angles{idx}.dist_zy   = rx_zy_dist{idx_zy};
            angles{idx}.dist_zx   = rx_zx_dist{idx_zx};
            angles{idx}.overlap_zy = rx_zy_overlap{idx_zy};
            angles{idx}.overlap_zx = rx_zx_overlap{idx_zx};
            angles{idx}.weight    = weight{idx_zy,idx_zx};
            
            % Image Array of angles
            angles{idx}.rx_zy_ar = rx_rot_zy;
            angles{idx}.rx_zx_ar = rx_rot_zx;
            angles{idx}.tx_zy_ar = tx_rot_zy;
            angles{idx}.tx_zx_ar = tx_rot_zx;
            
        end
    end
else % st.no_overlap == 0
   % Create output cell array
    idx=0;
    for idx_zy = 1:no_TX_zy
        for idx_zx = 1:no_TX_zx
            idx = idx+1;
            angles{idx}.rx_zx     = rx_zx{idx_zx};
            angles{idx}.rx_zy     = rx_zy{idx_zy};
            angles{idx}.rx_zx_idx = rx_zx_idx{idx_zx};
            angles{idx}.rx_zy_idx = rx_zy_idx{idx_zy};
            angles{idx}.tx_zy     = tx_rot_zy(idx_zy);
            angles{idx}.tx_zx     = tx_rot_zx(idx_zx);
            angles{idx}.weight    = ones(no_lines_per_tx_zy,no_lines_per_tx_zx);

            
            % Image Array of angles
            angles{idx}.rx_zy_ar = rx_rot_zy;
            angles{idx}.rx_zx_ar = rx_rot_zx;
            angles{idx}.tx_zy_ar = tx_rot_zy;
            angles{idx}.tx_zx_ar = tx_rot_zx;
            
        end
    end
end





%% Working loop to determine weight
function [weight ] = working_loop (idx_tx_zy, idx_tx_zx, ...
                                   no_TX_zx, no_TX_zy, ...
                                   no_lines_per_tx_zx, no_lines_per_tx_zy,...
                                   no_overlap, dist_src_mtr) 

% $$$ fprintf('zx: %02i  zy: %02i\n', idx_tx_zy, idx_tx_zx);
active_sources = (1:9);
%% remove not used sources at edges
if idx_tx_zy == 1
    for src = [1 4 7]
        active_sources(active_sources==src) = [];
    end
end
if idx_tx_zy == no_TX_zy
    for src = [3 6 9]
        active_sources(active_sources==src) = [];
    end
end
if idx_tx_zx == 1
    for src = [1 2 3]
        active_sources(active_sources==src) = [];
    end
end
if idx_tx_zx == no_TX_zx
    for src = [7 8 9]
        active_sources(active_sources==src) = [];
    end
end


%% matrix of number of sources per line
no_sources_small = ones(no_lines_per_tx_zy, no_lines_per_tx_zx);
for idx_src=active_sources
    % test for source
    switch idx_src
      case 1
        no_sources_small(1:no_overlap, 1:no_overlap) = 1 + no_sources_small(1:no_overlap, 1:no_overlap);        
      case 2
        no_sources_small(:,1:no_overlap) = 1 + no_sources_small(:,1:no_overlap);
      case 3
        no_sources_small(end:-1:end-no_overlap+1, 1:no_overlap) = 1 + no_sources_small(end:-1:end-no_overlap+1, 1:no_overlap);
      case 4
        no_sources_small(1:no_overlap, :) = 1 + no_sources_small(1:no_overlap, :);
      case 6
        no_sources_small(end:-1:end-no_overlap+1,:) = 1 + no_sources_small(end:-1:end-no_overlap+1,:);
      case 7
        no_sources_small(1:no_overlap, end:-1:end-no_overlap+1) = 1 + no_sources_small(1:no_overlap, end:-1:end-no_overlap+1);
      case 8
        no_sources_small(:,end:-1:end-no_overlap+1) = 1 + no_sources_small(:,end:-1:end-no_overlap+1);
      case 9
        no_sources_small(end:-1:end-no_overlap+1, end:-1:end-no_overlap+1) = 1 + no_sources_small(end:-1:end-no_overlap+1, end:-1:end-no_overlap+1);
    end  
end


%% expand to overlapping lines
no_sources = zeros(no_lines_per_tx_zy+no_overlap*2, no_lines_per_tx_zx+no_overlap*2);
no_sources(no_overlap+1:end-no_overlap, no_overlap+1:end-no_overlap) = no_sources_small;
% top left corner
no_sources(no_overlap:-1:1,no_overlap:-1:1) = no_sources_small(1:no_overlap,1:no_overlap);
% top right corner
no_sources(no_overlap:-1:1,end-no_overlap+1:end) = no_sources_small(1:no_overlap,end:-1:end-no_overlap+1);
% bottom left corner
no_sources(end-no_overlap+1:end,no_overlap:-1:1) = no_sources_small(end:-1:end-no_overlap+1,1:no_overlap);
% bottom right corner
no_sources(end-no_overlap+1:end,end-no_overlap+1:end) = no_sources_small(end:-1:end-no_overlap+1,end:-1:end-no_overlap+1);
% top center
no_sources(no_overlap:-1:1,no_overlap+1:end-no_overlap) = no_sources_small(1:no_overlap,1:end);
% bottom center
no_sources(end-no_overlap+1:end,no_overlap+1:end-no_overlap) = no_sources_small(end:-1:end-no_overlap+1,1:end);
% left center
no_sources(no_overlap+1:end-no_overlap,1:no_overlap) = no_sources_small(1:end, no_overlap:-1:1);
% right center
no_sources(no_overlap+1:end-no_overlap,end-no_overlap+1:end) = no_sources_small(1:end, end:-1:end-no_overlap+1);


%% Sum of distances per line within group
[dimy dimx] = size(no_sources);
dist_tot    = zeros(dimy,dimx);
weight_tot  = zeros(dimy,dimx);

source_tmp = [];
for idx_zy = 1:dimy
    for idx_zx = 1:dimx
        active_sources_tmp = active_sources;
        for idx = 1:no_sources(idx_zy,idx_zx)
            %idxx = idx
            [val_min idx_min] = min(dist_src_mtr(idx_zy,idx_zx,active_sources_tmp)); % get closest source
% $$$             source_tmp(idx_zy,idx_zx,idx) = active_sources_tmp(idx_min);             % save closest source
            active_sources_tmp(idx_min) = [];                                        % remove found source
% $$$             dist_tot(idx_zx,idx_zy)     = dist_tot(idx_zx,idx_zy) + val_min; 
            weight_tot(idx_zy,idx_zx)   = weight_tot(idx_zy,idx_zx) + 1/val_min;
        end
    end
end


%% remove non-active lines
center_src_dist = dist_src_mtr(:,:,5);
if idx_tx_zy == 1
    center_src_dist = center_src_dist(no_overlap+1:end, :);
    weight_tot      = weight_tot(no_overlap+1:end, :);
end
if idx_tx_zy == no_TX_zy
    center_src_dist = center_src_dist(1:end-no_overlap, :);
    weight_tot      = weight_tot(1:end-no_overlap, :);
end
if idx_tx_zx == 1
    center_src_dist = center_src_dist(:,no_overlap+1:end);
    weight_tot = weight_tot(:,no_overlap+1:end);
end
if idx_tx_zx == no_TX_zx
    center_src_dist = center_src_dist(:,1:end-no_overlap);
    weight_tot = weight_tot(:,1:end-no_overlap);
end


%% Calculate the actual weight
weight = (1./center_src_dist)./weight_tot;
weight(isnan(weight)) = 1;

% $$$ if idx_tx_zy == 1 && idx_tx_zx == 1
% $$$     figure;imagesc(weight_tot);title('wtot')
% $$$     figure;imagesc(center_src_dist);title('src_dist')
% $$$     figure;imagesc(weight);title('w')
% $$$ end

% $$$ weight2 = weight + fliplr(weight) + flipud(weight)  + fliplr(flipud(weight))
% $$$ ws_int = size(source_tmp)
% $$$ figure;imagesc(no_sources);title('sources')
% $$$ figure;imagesc(weight);title('weight')
% $$$ figure;imagesc(source_tmp(:,:,1));title('sources idx 1')
% $$$ figure;imagesc(source_tmp(:,:,2));title('sources idx 2')
% $$$ figure;imagesc(source_tmp(:,:,3));title('sources idx 3')
% $$$ figure;imagesc(source_tmp(:,:,4));title('sources idx 4')
% $$$ figure;imagesc(weight2);title('weight2')

% $$$ source_ar{idx_tx_zx, idx_tx_zy}     = source_tmp;
% $$$ no_sources_ar{idx_tx_zx, idx_tx_zy} = no_sources;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mfr_exploso_lines_setup.m ends here
