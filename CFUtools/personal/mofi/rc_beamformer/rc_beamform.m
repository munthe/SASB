
function bfm_data = rc_beamform(data, src_pos, drn_pos, points, c ,time_axis)

p = points;
num_drn = size(drn_pos,1);
t_start = time_axis(1);

%% beamforming
num_points = size(p,1);
d_src = zeros(num_points,1);
d_drn = zeros(num_points,num_drn);


disp 'calculating ToF'
A_src  = src_pos(1,:)';
B_src  = src_pos(2,:)';
AB_src = B_src-A_src;
t2= tic;
% TODO: Speed optimization
parfor idx_p = 1:num_points
    % from source to point
    P = p(idx_p,:)';
    AP_src = P-A_src;
    BP_src = P-B_src;
    s_src = dot(AP_src,AB_src)./norm(AB_src).^2;

    if s_src < 0
        %disp ([mat2str(idx_p) ' s_tx <0']);
        d_src(idx_p) = norm(AP_src);
    elseif s_src > 1
        %disp ([mat2str(idx_p) ' s_tx >0']);
        d_src(idx_p) = norm(BP_src);
    else
        d_src(idx_p) = norm(cross(AB_src, AP_src)) / norm(AB_src);
    end



    % from point to drain
    tmp = zeros(1,num_drn);
    for idx_drn = 1:num_drn
        A_drn = squeeze(drn_pos(idx_drn, 1,:));
        B_drn = squeeze(drn_pos(idx_drn, 2,:));
        AB_drn = B_drn-A_drn;
        AP_drn = P-A_drn;
        BP_drn = P-B_drn;
        
        s_drn = dot(AP_drn,AB_drn)./norm(AB_drn).^2;
        
        if s_drn < 0
            tmp(idx_drn) = norm(AP_drn);
        elseif s_drn > 1
            tmp(idx_drn) = norm(BP_drn);
        else
            tmp(idx_drn) = norm(cross(AB_drn, AP_drn)) / norm(AB_drn);
        end
    end
    d_drn(idx_p,:) = tmp;
end
toc(t2)
ToF = (repmat(d_src,[1 num_drn]) + d_drn)/c + t_start;




%% Interpolate data
disp 'Interpolating Data'
%data = linspace(0,400, 1000)'* ones(1,num_drn);
%data = t.rf_column;
%time_axis = t.p.time_axis(1:end-1);

%bfm_data = zeros(num_points,1);
%tmp = zeros(1,num_points);
line_rf = zeros(num_points,num_drn);
t1 = tic;
parfor idx_drn = 1:num_drn
    tmp = zeros(1,num_points);
    for idx_p = 1:num_points
        tmp(idx_p) = interp1(time_axis, data(:,idx_drn), ToF(idx_p,idx_drn));
    end
    line_rf(:,idx_drn) = tmp;
end
toc(t1)



%% Apodization
disp 'Apodizing Element Data'
apo = hanning(num_drn)';
apo_tot = repmat(apo, [num_points 1]);
line_rf_apo = line_rf .* apo_tot;



%% Summation
disp 'Summing Data'
bfm_data = squeeze(sum(line_rf_apo,2));
%bfm_data = squeeze(reshape(bfm_data, size(X)));





return

%% Tests 
test_idxs =0;
if test_idxs
    apo = zeros(num_x*num_y, 1);
    apo_tx = apo;
    apo_rx = apo;
    apo_tx(idx_tx) = 1;
    apo_rx(idx_rx) = 1;

    xdc_tx = mfr_emission('c',1540, 'type','dense_custom',...
                          'no_elm_x',num_x, ...
                          'no_elm_y',num_y, ...
                          'apo',     apo_tx,...
                          'no_act_elm',num_x*num_y);

    xdc_rx = mfr_emission('c',1540, 'type','dense_custom',...
                          'no_elm_x',num_x, ...
                          'no_elm_y',num_y, ...
                          'apo',     apo_rx,...
                          'no_act_elm',num_x*num_y);
    figure;xdc_tx.plot_apo;title('tx')
    figure;xdc_rx.plot_apo;title('rx')
end
