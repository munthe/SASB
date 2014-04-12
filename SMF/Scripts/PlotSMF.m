function [] = PlotSMF()

load SMFimage3000x192.mat image useCaseParams
RFdata = image;

Data_env = abs(RFdata/max(RFdata(:)));
% Data_env = abs(hilbert(RFdata));

[TGC,point,val_at_point] = CalculateTGC(Data_env,[]);
Data_env = Data_env./repmat(TGC(:),1,size(Data_env,2));

Data_lg = logcompress(Data_env);
Data_lg(Data_lg<-60) = -60;

figure(1);
plot(Data_lg(:,96))
figure(2);
imagesc(Data_lg)
colormap(gray)

end

function [TGC,point,val_at_point] = CalculateTGC(Data_env,nr_unique_points)
filter_length = 40;
point = [];
val_at_point = [];
% check if we are calculating TGC on a speckle image or not
% the difference between a speckle image and an image from a water phantom 
% is the range from the largest value to the mean value

max_data = max(Data_env(:));
mean_data = mean(Data_env(:));

ratio = max_data/mean_data

if(ratio < 200) % we are running TGC on a speckle image
    fprintf('TGC for speckle image is applied\n')
    % search for lower and upper bound
    temp = log10(Data_env(:));
    temp = temp(~isinf(temp));
    X = linspace(min(temp),max(temp),512);
    hist_sasb = hist(temp,X)./numel(Data_env);
    hist_sasb_cum = cumsum(hist_sasb);
    % lower boundary is chosen to be at 0.1
    filt_low = 10^X(find(hist_sasb_cum > 0.1,1,'first'));
    % higher boundary is chosen to be at 0.9
    filt_high = 10^X(find(hist_sasb_cum < 0.9,1,'last'));

    % 
    mask = zeros(size(Data_env));
    mask(Data_env < filt_high & Data_env > filt_low) = 1;

    Data_env(mask == 0) = nan;
    y = nanmedian(Data_env,2);
    x = 1:length(y);
    x = x(~isnan(y));
    y = y(~isnan(y));

    Y = interp1(x,y,1:size(Data_env,1),'spline',y(end))'; 
    Y = filter(ones(1,filter_length)/filter_length,1,Y);
    Y(1:filter_length) = Y(filter_length+1);
    TGC = Y./max(Y);
else
    fprintf('TGC for water phantom image is applied\n')
    
    filter_order = 30;
    nr_points = 600;
% extract data from the center of the image +- 20 columns
col_nr = 20;
temp = Data_env(:,floor(size(Data_env,2)/2)-col_nr:ceil(size(Data_env,2)/2+col_nr));
% temp = Data_env(:,78:92);

% remove any inf data points
temp(isinf(temp)) = min(temp(~isinf(temp)));
% max of the rows
t = max(temp,[],2);
% find the peaks through depth - noise can disturbe so we filter
t_filt = filter(ones(1,filter_order),1,t);
% set first 200 samples to min(t_filt) to remove noise
t_filt(1:200) = min(t_filt);
% detrend
p = polyfit(1:length(t_filt),t_filt',1);
x = polyval(p,1:length(t_filt));
t_filt = t_filt-x';
% now find the X strongest points 
[val,index] = sort(t_filt,'descend');
% now find the unique points
index = index(1:nr_points);

% now we make a binary vector to mask strong points from weak points
zer = zeros(1,length(t));
zer(index) = 1;

% now we want to close small gaps in the binary vector
for k = 1:length(zer)-1
   if(sum(zer(k+1:min(length(zer),k+30))) > 1)
       index = find(zer(k+1:min(length(zer),k+30)) == 1,1,'last');
       zer(k:k+index) = 1;
   end
end

% now we want to find how many unique strong points we have
bw = bwlabel(zer);

point = zeros(max(bw),1);
val_at_point = zeros(max(bw),1);

for k = 1:max(bw)
    % extract data for this one unique sector
    index_start = find(bw == k,1)-20;
    index_end = find(bw == k,1,'last')+20;
%     data = t_filt(max(index_start,1):min(index_end,length(t_filt)),:);
    data = t(max(index_start,1):min(index_end,length(t_filt)),:);
    % find maximum
    [val, index] = max(max(data,[],2));
    point(k) = index+index_start-1;
    val_at_point(k) = val;
end

% now only select the number of points the user want
[val_at_point index_sort] = sort(val_at_point,'descend');
point = point(index_sort);

val_at_point = val_at_point(point > 1);
point = point(point > 1);

if(~isempty(nr_unique_points))
    point = point(1:min(nr_unique_points,length(point)));
    val_at_point = val_at_point(1:min(nr_unique_points,length(point)));
    
    % now find the RF value
    for k = 1:min(nr_unique_points,length(point))
        % extract data for this one unique sector
        index_start = point(k)-40;
        index_end = point(k)+20;
%         data = temp(index_start:index_end,:);
        data = temp(index_start:index_end,:);
        % find maximum
        [val, index] = max(max(data,[],2));
        point(k) = index+index_start-1;
        val_at_point(k) = val;
%         figure,imagesc(data)
    end

end

[point index_sort] = sort(point,'ascend');
val_at_point = val_at_point(index_sort);


% x = 1:length(t);

% 
% figure(21)
% set(gcf,'position',[912   703   560   420])
% subplot(3,1,1)
% plot(t)
% hold on
% plot(x(point),t(point),'ro')
% hold off
% 
% subplot(3,1,2)
% plot(bw)
% 
% subplot(3,1,3)
% plot(zer)

% generate curve
if(length(point)~= 1)
    x = point; 
    y = val_at_point; 
    xi = 1:length(t); 
    yi = interp1(x,y,xi,'cubic');
%     yi = interp1(x,y,xi,'spline');
    % make sure no amplification
    yi(1:x(1)) = yi(find(xi >= x(1),1,'first'));
    yi(x(end):end) = yi(find(xi <= x(end),1,'last'));
else
    yi = repmat(val_at_point,1,size(Data_env,1));
end

TGC = yi(:)./max(yi);
val_at_point = val_at_point(:)./max(yi);
% val_at_point = val_at_point(:)./TGC(point);

% TGC = repmat(yi',1,size(Data_env,2));

end
end    