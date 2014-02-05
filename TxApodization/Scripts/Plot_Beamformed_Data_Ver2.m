function obj = Plot_Beamformed_Data_Ver2(varargin)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Check if number of arguments is correct                 %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(rem(nargin,2) > 0)
    error('Wrong number of arguments')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              Parse arguments                          %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figNr = [];
loadpath = [];
show_tgc = 0;
stageOfBF = [];
typeOfBF = [];
data_type = [];
db_range = 90;
graylevel = 8;
TGC = [];
RFscaling = [];
gain = [];
compression_scheme ='BK_muLaw_RTSC_VP_DRC_SEL_0';% 'Linear';
% compression_scheme = 'Linear';
nr_unique_points = [];
graylevel_eq = [];
filter_coeff = [];
windowtissueq = [];

RFdata = [];
useCaseParams = [];

for k = 1:2:nargin
    switch(lower(varargin{k}))
        case{'rfdata'}
             RFdata = varargin{k+1};
        case{'usecaseparams'}
            useCaseParams = varargin{k+1};
        case{'windowtissueq'}
            windowtissueq = varargin{k+1};
        case{'filter'}
            filter_coeff = varargin{k+1};
        case{'figure nr' , 'fig nr' , 'fig_nr' , 'figure_nr'}
            figNr = varargin{k+1};
        case{'loadpath'}
            loadpath = varargin{k+1};
        case{'show tgc'}
            if(strcmpi('show tgc',varargin{k+1}))
                show_tgc = 1;
            end
        case{'normalize'}
            normalize = varargin{k+1};
        case{'nr_unique_points'}
            nr_unique_points = varargin{k+1};
        case{'stageofbf' , 'stage'}
            stageOfBF = lower(varargin{k+1});
        case{'typeofbf' , 'type of beamformation'}
            typeOfBF = lower(varargin{k+1});
        case{'data type'}
            switch(lower(varargin{k+1}))
                case{'simulation','simulated','sim'}
                    data_type = 'sim';    
                case{'measurement','measured'}
                    data_type = 'profocus';    
                otherwise
                    error('Data type is unknown')
            end
        case{'db_range' , 'db' , 'dynamic range'}
            db_range = varargin{k+1};
        case{'graylevel'}
            graylevel = varargin{k+1};
        case{'tgc'}
            TGC = varargin{k+1};    
        case{'rf scaling' , 'rf_scaling' , 'rf data scaling','RFscaling','rfscaling'}
            RFscaling = varargin{k+1}; 
        case{'gain'}
            gain = varargin{k+1}; 
        case{'graylevel_eq'}
            graylevel_eq = varargin{k+1}; 
        case{'compression scheme' , 'compression'}
            compression_scheme = varargin{k+1};
        otherwise
            fprintf('This option is not available: %s\n',varargin{k});      
    end
end

if(isempty(stageOfBF))
    error('Beamformation stage was not specified')
end
if(isempty(typeOfBF))
    error('Type of beamformation was not specified')
end
if(isempty(data_type))
     error('Data type was not specified')
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              Load data                                %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(isempty(RFdata))
    switch(stageOfBF)
        case{'first','first stage'}

            try
                str = sprintf('%sFirst_Stage_RF_Data_Beamformed%s%s%s%sFirst_Stage_RF_Data_Beamformed.mat',loadpath,filesep,typeOfBF,data_type,filesep);   
                load(str,'RFdata','useCaseParams')
            catch
                load(loadpath,'RFdata','useCaseParams')
            end

        case{'second','second stage'}
            try
                str = sprintf('%sSecond_Stage_RF_Data_Beamformed%s%s%s%sSecond_Stage_RF_Data_Beamformed.mat',loadpath,filesep,typeOfBF,data_type,filesep); 
                load(str,'RFdata','useCaseParams','xmt')
            catch
                str = sprintf('%sSecond_Stage_RF_Data_Beamformeds%s%s%sSecond_Stage_RF_Data_Beamformed_1.mat',loadpath,filesep,typeOfBF,data_type,filesep); 
                load(loadpath,'RFdata','useCaseParams','xmt')
            end
        otherwise
            error('unknown stage of beamformation')
    end
else
    if(isempty(useCaseParams))
        error('missing usecaseparams file')        
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       Parse view                                      %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~isempty(windowtissueq))
    useCaseParams.scanparams(1).windowtissueq = windowtissueq;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       Calculate scan setup                            %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
line = CalculateScanSetup(RFdata,useCaseParams);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          RF data scaling                              %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RFscaling(isinf(RFscaling)) = 1;
RFscaling (RFscaling > 1) = 1;
if(~isempty(RFscaling))
    RFdata = RFdata.*RFscaling(1:size(RFdata,1),:);
end
RFdata(isnan(RFdata)) = 0;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                         Envelope detection                            %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for km = 1:1:405
useCaseParams.scanparams(1).scantype  = 10;
if(isreal(RFdata))
    if(useCaseParams.scanparams(1).scantype == 1)
        Data_env = RFdata;
    else
        Data_env = abs(hilbert(RFdata));
    end
else
    switch(useCaseParams.scanparams(1).scantype)
        case(1)
%             imgxz = zeros(size(RFdata,1),sqrt(size(RFdata,2)));
%             imgyz = zeros(size(RFdata,1),sqrt(size(RFdata,2)));
            for k = 1:size(RFdata,1) 
                d = reshape(abs(RFdata(k,:)),sqrt(size(RFdata,2)),sqrt(size(RFdata,2)));
                Data_c_scan_sum(k,:,:) = d;
%                 imgxz(k,:) = d(40,:);
%                 imgyz(k,:) = d(:,40);
            end

%             % xz
%             Data_env = abs(imgxz);
%             % xy
%             Data_env = abs(imgyz);

            % C scan
%             RFdata(1:200,:) = 0; % hack to remove scatter
            [val,id] = max(abs(RFdata),[],2);
            [dummy, depth_id] = max(val);
%             depth_id = 190;
%             depth_id = 811;

            
            
            Data_env = reshape(abs(RFdata(depth_id,:)),sqrt(size(RFdata,2)),sqrt(size(RFdata,2)));
            Data_c_scan_max = Data_env;
            Data_c_scan_sum = squeeze(sum(Data_c_scan_sum,1));
        case(4)
            % C scan
            Data_env = reshape(abs(RFdata(km,:)),sqrt(size(RFdata,2)),sqrt(size(RFdata,2)));
        otherwise
            Data_env = abs(RFdata);
    end
    
end
% save('SASB_c_scan_angl_2','Data_c_scan_sum','Data_c_scan_max')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                         Apply digital TGC                             %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~isempty(TGC)) % if not empty
    Data_env = Data_env./repmat(TGC(:),1,size(Data_env,2));
else
    if(useCaseParams.scanparams(1).scantype == 1)
    else
        [TGC,point,val_at_point] = CalculateTGC(Data_env,nr_unique_points);
        Data_env = Data_env./repmat(TGC(:),1,size(Data_env,2));
        if(~isempty(point))
            Data_env(Data_env > max(Data_env(point(1),:))) = max(Data_env(point(1),:));
        end
    end
end
obj.TGC = TGC;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                             Apply gain                                %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(isempty(gain))
    max_ = max(Data_env(:));
    min_ = min(Data_env(:));
    mask = Data_env < max_*0.8 & Data_env > min_ *1.2;
    med_ = median(Data_env(mask));
    gain = 10^(-20/10);

%     Data_env = Data_env./med_*gain;
    Data_env = Data_env./max(Data_env(:));
else
    max_ = max(Data_env(:));
    min_ = min(Data_env(:));
    mask = Data_env < max_*0.8 & Data_env > min_ *1.2;
    med_ = median(Data_env(mask));
    gain = 10^(gain/10);
    Data_env = Data_env./med_*gain;
end

disp(['max: ' num2str(max(Data_env(mask))) ' - median: ' num2str(median(Data_env(mask))) ' - mean: ' num2str(mean(Data_env(mask)))])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                            Log compress                               %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Data_lg,Data_gray] = logCompress(Data_env,compression_scheme, graylevel, db_range);

if(strcmp(compression_scheme,'BK_muLaw_RTSC_VP_DRC_SEL_0'))
    Data_lg = Data_lg - max(Data_lg(:));
end



switch([lower(typeOfBF) data_type])
    case{'drfprofocus'}
        obj.Tx_focus = useCaseParams.bfxmitparams(9).xmitfocus;
        obj.Tx_fnum = useCaseParams.bfxmitparams(9).xmitfnum;
        obj.Rx_focus = useCaseParams.bfrcvparams(3).rcvfocus;
        obj.Rx_fnum = useCaseParams.bfrcvparams(3).rcvfnum;
        obj.NrLines = useCaseParams.scanparams(2).stoplinenumq - useCaseParams.scanparams(2).startlinenumq + 1;
        obj.Method{1} = 'DRF';
        obj.Method{2} = 'measurement';
    case{'sasbprofocus'}
        obj.Tx_focus = useCaseParams.bfxmitparams(1).xmitfocus;
        obj.Tx_fnum = useCaseParams.bfxmitparams(1).xmitfnum;
        obj.Rx_focus = useCaseParams.bfrcvparams(1).rcvfocus;
        obj.Rx_fnum  = useCaseParams.bfrcvparams(1).rcvfnum;
        obj.NrLines = useCaseParams.scanparams(1).stoplinenumq - useCaseParams.scanparams(1).startlinenumq + 1;
        obj.Method{1} = 'SASB';
        obj.Method{2} = 'measurement';
        
    case{'drfsim'}
        obj.Tx_focus = useCaseParams.bfxmitparams(1).xmitfocus;
        obj.Tx_fnum = useCaseParams.bfxmitparams(1).xmitfnum;
        obj.Rx_focus = useCaseParams.bfrcvparams(1).rcvfocus;
        obj.Rx_fnum  = useCaseParams.bfrcvparams(1).rcvfnum;
        obj.NrLines = useCaseParams.scanparams(1).stoplinenumq - useCaseParams.scanparams(1).startlinenumq + 1;
        obj.Method{1} = 'DRF';   
        obj.Method{2} = 'simulation';   
    otherwise
%         obj.Tx_focus = useCaseParams.bfxmitparams(1).xmitfocus;
%         obj.Tx_fnum = useCaseParams.bfxmitparams(1).xmitfnum;
%        obj.Rx_focus = useCaseParams.bfrcvparams(1).rcvfocus;
        obj.Rx_fnum  = useCaseParams.bfrcvparams(1).rcvfnum;
        obj.NrLines = useCaseParams.scanparams(1).stoplinenumq - useCaseParams.scanparams(1).startlinenumq + 1;
%         obj.NrVS = xmt.no_lines;
        obj.Method{1} = 'SASB';
        obj.Method{2} = 'simulation';
        obj.Gain = gain;
end

try
% [frame, WindowTissueX , WindowTissueY] = ScanConvert(((Data_lg./db_range)+1).*2^16,...           % Data
%     useCaseParams.scanparams(1).scanareadscr.startlineorigin.x,...           % StartLineX
%     useCaseParams.scanparams(1).scanareadscr.startlineorigin.y,...           % StartLineY
%     useCaseParams.scanparams(1).scanareadscr.startlineangle,...              % StartLineAngle
%     useCaseParams.scanparams(1).startdepthq,...                  % StartDepth 
%     useCaseParams.scanparams(1).scanareadscr.stoplineorigin.x,...            % StopLineX
%     useCaseParams.scanparams(1).scanareadscr.stoplineorigin.y,...            % StopLineY
%     useCaseParams.scanparams(1).scanareadscr.stoplineangle,...               % StopLineAngle
%     useCaseParams.scanparams(1).startdepthq+line.length_lines,...% StopDepth
%     -.06,...                   % WindowTissueXMin 
%     useCaseParams.scanparams(1).windowtissueq.y_tismin,...                   % WindowTissueYMin 
%     0.06,...                   % WindowTissueXMax
%     useCaseParams.scanparams(1).windowtissueq.y_tismax);                     % WindowTissueYMax


[frame, WindowTissueX , WindowTissueY] = ScanConvert(((Data_lg./db_range)+1).*2^16,...           % Data
    useCaseParams.scanparams(1).scanareadscr.startlineorigin.x,...           % StartLineX
    useCaseParams.scanparams(1).scanareadscr.startlineorigin.y,...           % StartLineY
    useCaseParams.scanparams(1).scanareadscr.startlineangle,...              % StartLineAngle
    useCaseParams.scanparams(1).startdepthq,...                  % StartDepth 
    useCaseParams.scanparams(1).scanareadscr.stoplineorigin.x,...            % StopLineX
    useCaseParams.scanparams(1).scanareadscr.stoplineorigin.y,...            % StopLineY
    useCaseParams.scanparams(1).scanareadscr.stoplineangle,...               % StopLineAngle
    useCaseParams.scanparams(1).startdepthq+line.length_lines,...% StopDepth
    useCaseParams.scanparams(1).windowtissueq.x_tismin,...                   % WindowTissueXMin 
    useCaseParams.scanparams(1).windowtissueq.y_tismin,...                   % WindowTissueYMin 
    useCaseParams.scanparams(1).windowtissueq.x_tismax,...                   % WindowTissueXMax
    useCaseParams.scanparams(1).windowtissueq.y_tismax);                     % WindowTissueYMax


catch
    
    
[frame, WindowTissueX , WindowTissueY] = ScanConvert(((Data_lg./db_range)+1).*2^16,...           % Data
    useCaseParams.scanparams(1).scanareadscr.startlineorigin.y,...           % StartLineX
    useCaseParams.scanparams(1).scanareadscr.startlineorigin.x,...           % StartLineY
    useCaseParams.scanparams(1).scanareadscr.startlineangle,...              % StartLineAngle
    useCaseParams.scanparams(1).startdepthq,...                  % StartDepth 
    useCaseParams.scanparams(1).scanareadscr.stoplineorigin.y,...            % StopLineX
    useCaseParams.scanparams(1).scanareadscr.stoplineorigin.x,...            % StopLineY
    useCaseParams.scanparams(1).scanareadscr.stoplineangle,...               % StopLineAngle
    useCaseParams.scanparams(1).startdepthq+line.length_lines,...% StopDepth
    useCaseParams.scanparams(1).windowtissueq.x_tismin,...                   % WindowTissueXMin 
    useCaseParams.scanparams(1).windowtissueq.y_tismin,...                   % WindowTissueYMin 
    useCaseParams.scanparams(1).windowtissueq.x_tismax,...                   % WindowTissueXMax
    useCaseParams.scanparams(1).windowtissueq.y_tismax);                     % WindowTissueYMax
end

frame = ((frame./2^16)-1).*db_range;
if(~isempty(filter_coeff))
    frame = medfilt2(frame,filter_coeff);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                             Apply graylevel equalization              %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(isempty(graylevel_eq))
else
    v = linspace(-60,0,256);
    h = hist(frame(:),v);
    h = h./sum(h(:));
    h(1) = 0;
    graylevel_eq(1) = 0;
    [x,lag] = xcorr(graylevel_eq,h);

    [val_, index_]  = max(x);
    lag(index_)
    
    frame(frame(:) > v(10) & frame(:) < v(end)) = frame(frame(:) > v(10) & frame(:) < v(end)) +  lag(index_);
    
    h1 = hist(frame(:),v);
    h1 = h1./sum(h1(:));
% 
%     figure(199)
%     plot(h)
%     hold on
%     plot(h1,'r')
%     plot(graylevel_eq,'k')
%     hold off
%     
    [x,lag] = xcorr(graylevel_eq(10:end),h1(10:end));

    [val_, index_]  = max(x);
    
    lag(index_)
    
%     figure(4)
% plot(graylevel_eq(10:end))
% hold on
% plot(h(10:end),'r')
% hold off
end





frame(frame < -db_range) = -db_range;
% Visualize
if(isempty(figNr))
    f = figure;
else
    f = figure(figNr);
    clf(f)
end
%%

set(gcf,'PaperPositionMode','auto')
if(show_tgc)
    RECT_figure = [10 50 size(frame,2) + 50+285 size(frame,1)+160];
    RECT_axes1 = [80 80 size(frame,2) size(frame,1)];
    RECT_axes2 = [80+size(frame,2)+30 80 200 size(frame,1)];
else
    RECT_figure = [10 50 size(frame,2)+50+20 size(frame,1)*2+80];
    RECT_axes1 = [80 70 size(frame,2) size(frame,1)];
    
end

% set(f,'Position',RECT_figure)
axes1 = axes('YDir','reverse',...
    'units','pixel',...
    'Position',RECT_axes1,'FontName','Times New Roman');
if(strcmp(compression_scheme,'Linear'))
    db = linspace(-db_range,0,256);
%     imagesc(WindowTissueX*1000,WindowTissueY*1000,db(frame+1),'parent',axes1)
    
    if(useCaseParams.scanparams(1).scantype == 1)
%         imagesc(WindowTissueX*1000,WindowTissueY*1000,frame,'parent',axes1)
        imagesc((WindowTissueX-(pi/2))/(2*pi)*360,(WindowTissueY-pi/2)/(2*pi)*360,frame,'parent',axes1)
        set(axes1,'Position',[80 50 400 400])
        set(gcf,'Position',[1043         485         509         495])
        caxis([-50 0])

        return
    else
        imagesc(WindowTissueX*1000,WindowTissueY*1000,frame,'parent',axes1)
    end
    
    
else
    if(useCaseParams.scanparams(1).scantype == 1)
        imagesc((WindowTissueX-(pi/2))/(2*pi)*360,(WindowTissueY-pi/2)/(2*pi)*360,frame,'parent',axes1)
        set(axes1,'Position',[80 50 400 400])
        set(gcf,'Position',[1043         485         509         495])
        colormap(gray(2^graylevel))
        drawnow

        return
    else
        imagesc(WindowTissueX*1000,WindowTissueY*1000,frame,'parent',axes1)
    end
end
colormap(gray(2^graylevel))

if(show_tgc)
    axes2 = axes('YDir','reverse',...
    'units','pixel',...
    'Position',RECT_axes2,'FontName','Times New Roman');



    z = linspace(useCaseParams.scanparams(1).startdepthq,...
        useCaseParams.scanparams(1).startdepthq+line.length_lines,...
        length(TGC));

    plot(-20*log10(TGC),z*1000,'parent',axes2)
    hold on
    if(exist('val_at_point','var'))
        plot(-20*log10(val_at_point),z(point)*1000,'ro','parent',axes2)
    end
    axis(axes2,[-30 30 WindowTissueY(1)*1000 WindowTissueY(end)*1000])
    title('TGC','fontsize',12,'FontName','Times New Roman')
    xlabel('dB','fontsize',12,'FontName','Times New Roman')

    set(axes2,'YDir','reverse')
end
xlabel('Lateral position [mm]','fontsize',18,'parent',axes1)
ylabel('Axial position [mm]','fontsize',18,'parent',axes1)
set(gca,'fontsize',18)
axes(axes1)

% set(gcf,'Position',[10 50 1.5*340*1.1279 590*1.5])
% set(gca,'Position',[80 60 500*0.5765*1.5 520*1.5])
if(exist('axes2','var'))
    set(gcf,'Position',[ 10    50   675   885])
    set(axes2,'Position',[500*0.765*1.5 60 100*0.5765*1.5 520*1.5])
    tgc.TGC = TGC;
    if(exist('val_at_point','var'))
        tgc.point = point;
        tgc.val_at_point = val_at_point;
        tgc.z = z;
        set(gcf,'UserData',tgc)
    end
else
    set(gcf,'Position',[10 50 1.5*340*1.1279 590*1.5])
end


%%
% 
% clims = [0 db_range];
% imagesc(WindowTissueX*1000,WindowTissueY*1000,frame,clims)
% axis image
% xlabel('Lateral position [mm]','fontsize',16)
% ylabel('Axial position [mm]','fontsize',16)
% switch(obj.Method{1})
%     case{'SASB'}
%         title({[obj.Method{1} ' - ' obj.Method{2}],['Tx / Rx focus: ' num2str(obj.Tx_focus) ...
%                 ' / ' num2str(obj.Rx_focus) ' - Tx / Rx F#: ' num2str(obj.Tx_fnum) ...
%                ' / ' num2str(obj.Rx_fnum)],['Nr. lines: ' num2str(obj.NrLines) ' - dynamic range ' num2str(db_range) 'dB'],...
%                },'fontsize',16);
%     case{'DRF'}
%         title({[obj.Method{1} ' - ' obj.Method{2}],['Tx focus: ' num2str(obj.Tx_focus) ...
%                 ' - Tx / Rx F#: ' num2str(obj.Tx_fnum) ...
%                ' / ' num2str(obj.Rx_fnum)],['Nr. lines: ' num2str(obj.NrLines) ' - dynamic range ' num2str(db_range) 'dB'],...
%                },'fontsize',16);    
% end
% 
% set(gcf,'PaperPositionMode','auto')
% 
% 
% 
% if(show_tgc == 1)
%     set(gca,'units','pixel')
%     set(f,'Position',[9          61         887        1061])
%     ax1_h = gca;
%     ax1 = get(ax1_h,'position');
%     ax = axes;
%     set(ax,'units','pixels');
%     set(ax,'position',round([ax1(1)+ax1(3)+50 ax1(2) 150 ax1(4)]))
% 
%     z = linspace(useCaseParams.scanparams(1).scanareadscr.startdepth,...
%         useCaseParams.scanparams(1).scanareadscr.startdepth+line.length_lines,...
%         length(t));
%     plot(-y,z(x),'o',-yi,z(xi))
%     axis(ax,[-30 30 WindowTissueY(1) WindowTissueY(end)])
%     axis(ax,'ij')
%     title('TGC','fontsize',12)
%     xlabel('dB','fontsize',12)
%     axes(ax1_h);
% end



%% Filename
switch(useCaseParams.bfrcvparams(1).rcvapodishape)
    case(0) %{'rect','boxcar'}
        obj.apotype = 'Rectwin';
    case(1) %{'hamming','hamm'}
        obj.apotype = 'Hamming';   
    case(2) %{'gauss'}
        obj.apotype = 'Gauss';    
    case(3) %{'hanning','hann'}
        obj.apotype = 'Hann';
    case(4) %{'black','blackman'}
        obj.apotype = 'Blackman';
    case(5) %{'bartlett'}
        obj.apotype = 'Bartlett';
end

switch(obj.Method{1})
    case{'SASB'}
%     obj.filename = [obj.Method{1} '_' obj.Method{2} '_VS_' num2str(obj.Tx_focus) '_VS_Fnum_' num2str(obj.Tx_fnum) ...
%         '_NrVS_' num2str(obj.NrVS) '_Apo_' obj.apotype '_Bmode'];
%     title(strrep(obj.filename,'_',' '))
    case{'DRF'}
    obj.filename = [obj.Method{1} '_' obj.Method{2} '_Tx_focus_' num2str(obj.Tx_focus) '_Fnum_' num2str(obj.Tx_fnum) '_Rx_Fnum_' num2str(obj.Rx_fnum)...
        '_NrLines_' num2str(obj.NrLines) '_Apo_' obj.apotype '_Image'];

end





% end
% obj.filename = strrep(obj.filename,'.','_');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                         Internal functions                            %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [line] = CalculateScanSetup(RFdata,useCaseParams)
if(~isempty(RFdata))
    % Calculate scan line reference positions
    [rcv.scanline_ref_point rcv.scanline_direction rcv.scanline_angle] = calc_scanline_position(useCaseParams);

    line.dr              = useCaseParams.scanparams(1).c_sound/useCaseParams.bfrcvparams(1).smpfreq/2;
    line.direction       = rcv.scanline_direction;

    % Calculate layer thickness
    aclayerthickness = useCaseParams.acmodparams(1).layerthickness1 + ...
                       useCaseParams.acmodparams(1).layerthickness2 + ...
                       useCaseParams.acmodparams(1).layerthickness3;
    % Calculate line origin which is normaly at the transducer surface and not
    % on the element surface
    line.origin          = rcv.scanline_ref_point + line.direction*aclayerthickness;
    % In bft3 the line must start where the first sample is wanted.
    % So we correct the origin with startdepthq along the direction of the line
    line.origin          = line.origin +  line.direction*useCaseParams.scanparams(1).startdepthq;

    line.length_lines    = (size(RFdata,1)-1)*line.dr;                
    line.N_points_in_line = size(RFdata,1);
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
    
    
    
    