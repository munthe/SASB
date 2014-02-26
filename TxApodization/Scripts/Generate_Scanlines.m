%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   useCaseParams       Usecase structure
%   no_active_elements  maximum number of active elements
%   symmetric           'symmetric'/'nonsymmetric' apodization
%   mode                'receive'/'transmit'
%   Optional:
%       * varargin{1}: debug option. 'debug' if plot is wanted
% OUTPUT
%   rcv                         structure
%    .no_lines                  number of scan lines 
%    .focus                     focus depth relative to scanline_ref_point
%    .fnum                      F#
%    .scanline_direction        vector of unitvector directions [lineNr,3]     
%    .scanline_angle            vectof of angles (angle)[lineNr,1]     
%    .scanline_ref_point        vector of reference point [lineNr,3]
%    .focus_point               vector of focus points [lineNr,3]
%    .dist_scanline_to_elements matrix of distances [lineNr,elementNr]
%    .valid                     matrix of valid elements [lineNr,elementNr]
%    .apo                       matrix of apo. values [lineNr,elementNr]
%    .apodishape                apodization shape number from the usecase
%    .delays                    time delay profile
%    .delay_offsets             time delay offset to use with Field II
%
% DEPENDENCIES
%   M-files: transducer class, calc_trm_delays, calc_scanline_position
%
% VERSION		
%   v1  2010-05-26
% AUTHOR    Martin Christian Hemmsen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rcv = Generate_Scanlines(useCaseParams,no_active_elements,symmetric,mode,varargin)

debug = 0;
if(nargin > 4)
    debug = varargin{1};
end
        

switch(mode)
    case{'receive','rcv'}
        mode = 'bfrcvparams';
        focus_mode = 'rcvfocus';
        fnum_mode = 'rcvfnum';
        apodshape_mode = 'rcvapodishape';
        apodigausswidth_mode = 'rcvapodigausswidth';
        apodilevels_mode = 'rcvapodilevels'
    case{'transmit','xmt'}
        mode = 'bfxmitparams';
        focus_mode = 'xmitfocus';
        fnum_mode = 'xmitfnum';    
        apodshape_mode = 'xmitapodishape';
        apodigausswidth_mode = 'xmitapodigausswidth';
        apodilevels_mode = 'xmitapodilevels'
end
% calculate nr of lines

switch(useCaseParams.scanparams(1).scantype)
    case(0) % 1d imaging
        % define transducer
        transducerParams = transducer('id','unspecified',...
                        'kerf',useCaseParams.acmodparams(1).pitch/10,...
                        'pitch',useCaseParams.acmodparams(1).pitch,...
                        'aclayerthickness',useCaseParams.acmodparams(1).layerthickness1 + ...
                                           useCaseParams.acmodparams(1).layerthickness2 + ...
                                           useCaseParams.acmodparams(1).layerthickness3,...
                        'ROC',useCaseParams.acmodparams(1).shellradius, ...
                        'nr_elements_x',useCaseParams.acmodparams(1).elements,...
                        'nr_elements_y',0);
         rcv.no_lines = (useCaseParams.scanparams(1).stoplinenumq - useCaseParams.scanparams(1).startlinenumq + 1);    
        % define transducer
%         transducerParams = transducer('id','unspecified',...
%                         'kerf',useCaseParams.acmodparams(1).pitch/10,...
%                         'pitch',useCaseParams.acmodparams(1).pitch,...
%                         'aclayerthickness',useCaseParams.acmodparams(1).layerthickness1 + ...
%                                            useCaseParams.acmodparams(1).layerthickness2 + ...
%                                            useCaseParams.acmodparams(1).layerthickness3,...
%                         'ROC',useCaseParams.acmodparams(1).shellradius, ...
%                         'nr_elements_x',useCaseParams.acmodparams(1).elements,...
%                         'nr_elements_y',useCaseParams.acmodparams(1).elements);
                
    case(1) % 2d imaging
        rcv.no_lines = (useCaseParams.scanparams(1).stoplinenumq - useCaseParams.scanparams(1).startlinenumq + 1)^2;    
        % define transducer
        transducerParams = transducer('id','unspecified',...
                        'kerf',useCaseParams.acmodparams(1).pitch/10,...
                        'pitch',useCaseParams.acmodparams(1).pitch,...
                        'aclayerthickness',useCaseParams.acmodparams(1).layerthickness1 + ...
                                           useCaseParams.acmodparams(1).layerthickness2 + ...
                                           useCaseParams.acmodparams(1).layerthickness3,...
                        'ROC',useCaseParams.acmodparams(1).shellradius, ...
                        'nr_elements_x',useCaseParams.acmodparams(1).elements,...
                        'nr_elements_y',useCaseParams.acmodparams(1).elements);
    case(2) % xz imaging
        rcv.no_lines = (useCaseParams.scanparams(1).stoplinenumq - useCaseParams.scanparams(1).startlinenumq + 1);    
        % define transducer
        transducerParams = transducer('id','unspecified',...
                        'kerf',useCaseParams.acmodparams(1).pitch/10,...
                        'pitch',useCaseParams.acmodparams(1).pitch,...
                        'aclayerthickness',useCaseParams.acmodparams(1).layerthickness1 + ...
                                           useCaseParams.acmodparams(1).layerthickness2 + ...
                                           useCaseParams.acmodparams(1).layerthickness3,...
                        'ROC',useCaseParams.acmodparams(1).shellradius, ...
                        'nr_elements_x',useCaseParams.acmodparams(1).elements,...
                        'nr_elements_y',useCaseParams.acmodparams(1).elements);
    case(3) % yz imaging
        rcv.no_lines = (useCaseParams.scanparams(1).stoplinenumq - useCaseParams.scanparams(1).startlinenumq + 1);    
        % define transducer
        transducerParams = transducer('id','unspecified',...
                        'kerf',useCaseParams.acmodparams(1).pitch/10,...
                        'pitch',useCaseParams.acmodparams(1).pitch,...
                        'aclayerthickness',useCaseParams.acmodparams(1).layerthickness1 + ...
                                           useCaseParams.acmodparams(1).layerthickness2 + ...
                                           useCaseParams.acmodparams(1).layerthickness3,...
                        'ROC',useCaseParams.acmodparams(1).shellradius, ...
                        'nr_elements_x',useCaseParams.acmodparams(1).elements,...
                        'nr_elements_y',useCaseParams.acmodparams(1).elements);
end

    



rcv.focus = useCaseParams.(mode)(1).(focus_mode);
rcv.fnum = useCaseParams.(mode)(1).(fnum_mode);


if(rcv.focus(1) == 0) % in this case DRF focusing is assumed and we set focus 
                   % to some distance for calculation of distance between 
                   % elements and scanlines
    rcv.focus = 1;
    dyn_focus = 'enabled';
else
    dyn_focus = 'disabled';
end

if(length(rcv.focus) ~= rcv.no_lines)
    rcv.focus = repmat(rcv.focus, rcv.no_lines,1);
    rcv.fnum = repmat(rcv.fnum, rcv.no_lines,1);
end

% based on F# and focus calculate desired width of aperture
Width_desired = abs((rcv.focus./rcv.fnum));
% round width_desired up to match an equal nr of elements
Width_desired = ceil(Width_desired./useCaseParams.acmodparams(1).pitch)*useCaseParams.acmodparams(1).pitch;


% Calculate scan line ref positions
[rcv.scanline_ref_point rcv.scanline_direction rcv.scanline_angle] = calc_scanline_position(useCaseParams);


switch(mode)
    case{'bfrcvparams','receive','rcv'}
        % Calculate receive sample positions
        rcv.sample_origin          = rcv.scanline_ref_point + rcv.scanline_direction*transducerParams.aclayerthickness;
        % In bft3 the line must start where the first sample is wanted.
        % So we correct the origin with startdepthq along the direction of the line
        rcv.sample_origin          = rcv.sample_origin +  rcv.scanline_direction*useCaseParams.scanparams(1).startdepthq;

        rcv.dr = useCaseParams.scanparams(1).c_sound/useCaseParams.bfrcvparams(1).smpfreq/2;

    case{'bfxmitparams','transmit','xmt'}
end


rcv.focus_point = rcv.scanline_ref_point + repmat(rcv.focus,1,3).*rcv.scanline_direction;



% figure(3)
%  plot3(rcv.scanline_ref_point(:,1),rcv.scanline_ref_point(:,2), rcv.scanline_ref_point(:,3),'.')
% hold on
% plot3(rcv.focus_point(:,1),rcv.focus_point(:,2), rcv.focus_point(:,3),'bx')
% 
% for k = 1:1:rcv.no_lines
%     plot3([rcv.scanline_ref_point(k,1) rcv.focus_point(k,1)],[rcv.scanline_ref_point(k,2) rcv.focus_point(k,2)],[rcv.scanline_ref_point(k,3) rcv.focus_point(k,3)],'r')
% end
% hold off
% 

% calculate distance from scan line to elements
x0 = transducerParams.element_positions;

% 
% Geo = transducerParams.ReadGeoDebugParameters;
% 
% elementpos = flipud([Geo.element_position_table(:,1) Geo.element_position_table(:,2)]);
% scanline = [Geo.scan_line_table(:,2) Geo.scan_line_table(:,3)];
% 
% figure,plot( transducerParams.element_positions(:,1), transducerParams.element_positions(:,3),'.')
% hold on
% plot(rcv.scanline_ref_point(:,1),rcv.scanline_ref_point(:,3),'or')
% legend('element position (calculated)','scanline ref point (calculated')
% 
% 
% figure,plot(elementpos(:,1),elementpos(:,2),'or')
% hold on
% plot( transducerParams.element_positions(:,1), transducerParams.element_positions(:,3),'.')
% l = legend('element position (GeoDebug)','element position (calculated');
% set(l,'fontsize',12)
% set(gca,'fontsize',12)
% xlabel('Lateral displacement [mm]','fontsize',14)
% ylabel('Axial displacement [mm]','fontsize',14)
% 
% err = sqrt((transducerParams.element_positions(:,1)- elementpos(:,1)).^2 + ...
%           (transducerParams.element_positions(:,3)- elementpos(:,2)).^2);
%            
% figure,plot(err,'-')
% xlabel('Element nr','fontsize',14)
% ylabel('Error','fontsize',14)
% axis tight
% 
% figure,plot( scanline(:,1), scanline(:,2),'.')
% hold on
% plot(rcv.scanline_ref_point(:,1),rcv.scanline_ref_point(:,3),'or')
% l = legend('scanline ref point (GeoDebug)','scanline ref point (calculated');
% set(l,'fontsize',12)
% set(gca,'fontsize',12)
% xlabel('Lateral displacement [mm]','fontsize',14)
% ylabel('Axial displacement [mm]','fontsize',14)
% 
% 
% err = sqrt((rcv.scanline_ref_point(:,1)- scanline(:,1)).^2 + ...
%                (rcv.scanline_ref_point(:,3)- scanline(:,2)).^2);
%            
% figure,plot(err,'-')
% xlabel('Lateral displacement [mm]','fontsize',14)
% ylabel('Axial displacement [mm]','fontsize',14)
% title('Error','fontsize',16)

% old version - do not delete
%  rcv.dist_scanline_to_elements = zeros(rcv.no_lines,useCaseParams.acmodparams(1).elements^2);
rcv.dist_scanline_to_elements = zeros(rcv.no_lines,useCaseParams.acmodparams(1).elements);
for k = 1:rcv.no_lines
    % http://www.mathworks.se/support/solutions/en/data/1-1BYSR/
    x1 = rcv.focus_point(k,:);
    x2 = rcv.scanline_ref_point(k,:);   
    a = x1-x2;    
    b = x0 - repmat(x2,size(x0,1),1);
    rcv.dist_scanline_to_elements(k,:) = sum(abs(cross(repmat(a,size(x0,1),1),b)).^2,2).^(1/2) / norm(a);
end


% find valid elements
if(useCaseParams.scanparams(1).scanareadscr.startlineorigin.y ~= 0 & ...
        useCaseParams.scanparams(1).scanareadscr.startlineorigin.x ~= 0) % we assume 2d imaging
    rcv.valid = zeros(size(rcv.dist_scanline_to_elements));
    for k = 1:rcv.no_lines
        [value,index] = sort((rcv.dist_scanline_to_elements(k,:)));

        % enforce symmetric aperture arround scanline
        if(strcmpi(symmetric,'symmetric'))
            % check if number of valid elements is larger then
            % no_active_elements-1
            temp = value(1:no_active_elements-1) <= Width_desired(k)/2+eps;
            rcv.valid(k,index(temp)) = 1;
        else
            % check if number of valid elements is larger then no_active_elements
            temp = value(1:no_active_elements) <= Width_desired(k)/2+eps;
            rcv.valid(k,index(temp)) = 1;
        end
    end  
    
else
    rcv.valid = zeros(rcv.no_lines,useCaseParams.acmodparams(1).elements);
    for k = 1:rcv.no_lines
        [value,index] = sort((rcv.dist_scanline_to_elements(k,:)));

        % enforce symmetric aperture arround scanline
        if(strcmpi(symmetric,'symmetric'))
            % check if number of valid elements is larger then
            % no_active_elements-1
            temp = value(1:no_active_elements-1) <= Width_desired(k)/2+eps;
            rcv.valid(k,index(temp)) = 1;
        else
            % check if number of valid elements is larger then no_active_elements
            temp = value(1:no_active_elements) <= Width_desired(k)/2+eps;
            rcv.valid(k,index(temp)) = 1;
        end
    end   
end
% define apodization function
rcv.apodishape = useCaseParams.(mode)(1).(apodshape_mode);
switch rcv.apodishape
    case {0} %'boxcar','ones','none','rect','rectangular','square'
        apod_fun = @(n,beta) (n <= 1);
        rcv.apod_fun_name = 'Rectwindow';
    case {1} % 'hamming','hamm'
        apod_fun = @(n,beta) (0.54-0.46*cos(2.*n*pi)).*(n <= 1);
        rcv.apod_fun_name = 'Hamming';
    case {2} % 'gauss'
        apod_fun = @(n,beta) exp(-0.5*((n-0.5)/(beta*0.5)).^2).*(n <= 1);
        rcv.apod_fun_name = 'Gaussian';
    case {3} % 'hanning','hann'
        apod_fun = @(n,beta) (0.5*(1-cos(2.*n*pi))).*(n <= 1);
        rcv.apod_fun_name = 'Hann';
    case {4} % 'black','blackman'
        apod_fun = @(n,beta) ((1-0.16)/2-1/2*cos(2.*n*pi)+0.16/2*cos(2*2.*n*pi)).*(n <= 1);
        rcv.apod_fun_name = 'Blackman';
    case {5} % 'bart','bartlett'
        apod_fun = @(n,beta) (0.62-0.48*abs(n-1/2)-0.38*cos(2.*n*pi)).*(n <= 1);
        rcv.apod_fun_name = 'Bartlett-Hann';
    case {6} % 'tukey' (1/2*(1-cos(pi*(2*(n-0.5)/(beta))))).*((n-0.5) <= beta/2) + ...
        apod_fun = @(n,beta) ...
                             1.*((n-0.5) <= 0.5*(1-beta/2)) + ...
                             (1/2*(1+cos(pi*(2*(n-0.5)/(beta*0.5)-2/beta+1)))).*(0.5*(1-beta/2) < (n-0.5)).*(n-0.5 <= 0.5);
        rcv.apod_fun_name = 'Tukey';
    otherwise
        error('Wrong apodization type');
end
rcv.apo = zeros(size(rcv.dist_scanline_to_elements));
for k = 1:rcv.no_lines
    n = rcv.dist_scanline_to_elements(k,rcv.valid(k,:)==1)./(Width_desired(k))+0.5;
    rcv.apo(k,rcv.valid(k,:)==1) = apod_fun(n,useCaseParams.(mode)(1).(apodigausswidth_mode));
end
% Quantization of the apodization
if length( useCaseParams.(mode)(1).(apodilevels_mode) ) > 1
    rcv.apo = quantization(rcv.apo,useCaseParams.(mode)(1).(apodilevels_mode));
end
     
switch(dyn_focus)
    case{'disabled'}
        % calculate time delays
        rcv.delays = zeros(size(rcv.dist_scanline_to_elements));
        element_positions = transducerParams.element_positions;
        
        for k = 1:rcv.no_lines
            if(sum(rcv.valid(k,:) == 1) > 0)
            [rcv.delays(k,rcv.valid(k,:) == 1) rcv.delay_offsets(k)] = calc_trm_delays(...
                element_positions(rcv.valid(k,:) == 1,:),...
                rcv.focus_point(k,:),...
                useCaseParams.scanparams(1).c_sound,...
                rcv.scanline_ref_point(k,:));
            end

        end
    case{'enabled'}
        rcv.delays = zeros(size(rcv.dist_scanline_to_elements));   
end
%%
if(debug == 1)
    figure
    for k = 1:rcv.no_lines % only works as a sweeping function for plotting
        if(k == 1) % draw all elements
            transducerParams.plot_elements;
            hold on
            focus_point_plot =  plot(rcv.focus_point(k,1)*1000,rcv.focus_point(k,3)*1000,'xb','MarkerSize',10,'LineWidth',2,'DisplayName','Focus point');
            hold off
            % plot valid elements
            transducerParams.plot_elements(find(rcv.valid(k,:) == 1),[1 0 0],'Valid elements')
            hold on
            ref_point_plot = plot(rcv.scanline_ref_point(k,1)*1000,rcv.scanline_ref_point(k,3)*1000,'xk','MarkerSize',10,'LineWidth',2,'DisplayName','Scanline ref. point');
            % plot scan lines
            scan_line_plot = plot([rcv.scanline_ref_point(k,1) rcv.focus_point(k,1)]*1000,[rcv.scanline_ref_point(k,3) rcv.focus_point(k,3)]*1000,'--','color',[0 0 0],'linewidth',2,'DisplayName','Active scanline');
            hold off
            l = legend;
            set(l,'Location','EastOutside')
            axis image
            axis([min(min(rcv.focus_point(:,1)),min(transducerParams.element_positions(:,1))) max(max(rcv.focus_point(:,1)),max(transducerParams.element_positions(:,1))) min(-0.005,min(rcv.scanline_ref_point(:,3))) max(0.015,max(rcv.focus(:)))]*1000)
        else       % color the elements that is not active blue
            transducerParams.plot_elements(find((rcv.valid(k,:) == 0 & rcv.valid(k-1,:) == 1) == 1),[0 0 1],'');
            % plot focus points
            set(focus_point_plot,'XData',rcv.focus_point(k,1)*1000)
            set(focus_point_plot,'YData',rcv.focus_point(k,3)*1000)
            % plot valid elements
            transducerParams.plot_elements(find((rcv.valid(k,:) == 1 & rcv.valid(k-1,:) == 0) == 1),[1 0 0],'');
            % plot ref point
            set(ref_point_plot,'XData',rcv.scanline_ref_point(k,1)*1000)
            set(ref_point_plot,'YData',rcv.scanline_ref_point(k,3)*1000)
            % plot scan lines
            set(scan_line_plot,'XData',[rcv.scanline_ref_point(k,1) rcv.focus_point(k,1)]*1000)
            set(scan_line_plot,'YData',[rcv.scanline_ref_point(k,3) rcv.focus_point(k,3)]*1000)
           
        end
        
        drawnow
        
    end 
end
