% Varargin:
%   * beamformation_method:
%        * Single Receive Focus:  'sasb first stage' or 'srf'
%        * SASB second stage:     'sasbprofocus' or 'sasbsim'
%        * Dynamic Recieve Focus: 'drf'
%   * apodization: 
%        * Scalar: 1 Rect, 2 Han, 3 Ham, 4 Black, 5 Bartlett
%        * Vector: Specified in ReadMidLevelApodizationProfiles
%   * Delay profile:
%        * Vector: Specified in ReadMidLevelApodizationProfiles
%   * usecase
%   * loadpath
%   * 
%   * 
%   * 
function BeamformationVer4(varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Check if number of arguments is correct                 %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(rem(nargin,2) > 0)
    error('Wrong number of arguments')
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              Parse arguments                          %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
type = [];
bft3.apod.type = 1; % default rect window
useCaseParams = [];
loadpath = [];
debug = 0;
display = 1;
no_active_elements = 64;
symmetric = 'symmetric';
scanlines = [];
dly_matrix = [];
phase_profile = [];
generate_scaling = false;
savepath = [];
data = [];
tstart = 0;
xmt = [];
for k = 1:2:nargin
    switch(lower(varargin{k}))
        case{'generate_scaling'}
            generate_scaling = lower(varargin{k+1});
        case{'data'}
            data = varargin{k+1};
        case{'tstart'}
            tstart = varargin{k+1};

        case{'xmt'}
            xmt = varargin{k+1};
        case{'TransducerDescription'}
            TransducerDescription = varargin{k+1};
        case{'beamformation_method'}
            type = lower(varargin{k+1});
        case{'savepath'}
            savepath = varargin{k+1};
        case{'apodization'}
            if(ischar(varargin{k+1}))
                switch(lower(varargin{k+1}))
                    case{'rect','boxcar'}
                        useCaseParams.bfrcvparams(1).rcvapodishape = 0;
                    case{'hamming','hamm'}
                        useCaseParams.bfrcvparams(1).rcvapodishape = 1;
                    case{'gauss'}
                        useCaseParams.bfrcvparams(1).rcvapodishape = 2;
                    case{'hanning','hann'}
                        useCaseParams.bfrcvparams(1).rcvapodishape = 3;
                    case{'black','blackman'}
                        useCaseParams.bfrcvparams(1).rcvapodishape = 4;
                    case{'bart','bartlett'}
                        useCaseParams.bfrcvparams(1).rcvapodishape = 5;    
                    otherwise
                        error('Undefined apodization type')
                end
            else
                bft3.apod.type = varargin{k+1};
            end
        case{'delay_compensation_profile'}
            phase_profile = varargin{k+1};
        case{'delay'}
            dly_matrix = varargin{k+1};
        case{'usecase'}
            useCaseParams = varargin{k+1};
        case{'loadpath'}
            loadpath = varargin{k+1};
        case{'debug'}
            if(ischar(varargin{k+1}))
                switch(lower(varargin{k+1}))
                    case{'off'}
                        debug = 0;
                    case{'on'}
                        debug = 1;
                    otherwise
                        error('Undefined debug type')
                end
            else
                debug = varargin{k+1};
            end
        case{'symmetric'}
             if(ischar(varargin{k+1}))
                switch(lower(varargin{k+1}))
                    case{'symmetric'}
                        symmetric = 'symmetric';
                    case{'nonsymmetric'}
                        symmetric = 'nonsymmetric';
                    otherwise
                        error('Undefined symmetri value')
                end
             end
        case{'scanlines'}
            scanlines = varargin{k+1};
        case{'display'}
            display = varargin{k+1};
        otherwise
            fprintf('This option is not available: %s\n',varargin{k});      
    end    
end

if(isempty(type))
    error('Type was not specified')
end

if(isempty(useCaseParams))
    error('UseCase was not specified')
end

if(isempty(loadpath))
    error('Loadpath was not specified')
end
if(isempty(savepath))
    savepath = loadpath;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                   Parse useCaseParams for apodization                 %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bft3.apod.matrix = []; % default this is empty. 
if(max(size(bft3.apod.type)) == 1)
    switch(useCaseParams.bfrcvparams(1).rcvapodishape)
        case(0) %{'rect','boxcar'}
            bft3.apod.type = 'Rectwin';
        case(1) %{'hamming','hamm'}
            bft3.apod.type = 'Hamming';   
        case(2) %{'gauss'} 
            bft3.apod.type = 'Gaussian';    
        case(3) %{'hanning','hann'}
            bft3.apod.type = 'Hann';
        case(4) %{'black','blackman'}
            bft3.apod.type = 'Blackman';
        case(5) %{'bartlett'}
            bft3.apod.type = 'Bartlett';
    end
else
    bft3.apod.matrix = bft3.apod.type;
    bft3.apod.type = 'user specified apodization';
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Switch on the different beamformation types                   %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch(lower(type))
    case {'srf'}
        % remove any old directory
        if(isdir([loadpath 'First_Stage_RF_Data_Beamformed' filesep 'sasbsim' ]))
           rmdir([loadpath 'First_Stage_RF_Data_Beamformed' filesep 'sasbsim'],'s')
        end
        mkdir([loadpath 'First_Stage_RF_Data_Beamformed' filesep 'sasbsim'])
        
        xmt = Generate_Scanlines(useCaseParams, no_active_elements, symmetric, 'transmit',debug); 
        
        rcv = Generate_Scanlines(useCaseParams, no_active_elements, symmetric, 'receive',debug); 
        
        if(~isempty(bft3.apod.matrix)) % if user specified matrix overload
                                       % Generate_Scanlines apodization
            rcv.apo = bft3.apod.matrix;
            rcv.apod_fun_name = 'user specified apodization';
            if(display == 1)
                fprintf('Generate_Scanlines apodization is overloaded\n')
            end
        end
        
        if(~isempty(dly_matrix)) % if user specified delay overload
            rcv.delays = dly_matrix;
            rcv.delay_offsets = max(dly_matrix,[],2)';
            rcv.delays_modified = 'user specified delays';
            if(display == 1)
                fprintf('Generate_Scanlines delays are overloaded\n')
            end
        end
        
        bft3.apod.f_number    = useCaseParams.bfrcvparams(1).rcvfnum;

    case {'sasbprofocus','sasbsim','sasbsim_new'}   
        % remove any old directory
        if(isdir([loadpath 'Second_Stage_RF_Data_Beamformed' filesep lower(type) ]))
           rmdir([loadpath 'Second_Stage_RF_Data_Beamformed' filesep lower(type)],'s')
        end
        mkdir([loadpath 'Second_Stage_RF_Data_Beamformed' filesep lower(type)])
        
        if(isempty(xmt))
            % Load parameters from the 1st stage beamformer
            str = sprintf('%sFirst_Stage_RF_Data_Beamformed%ssasbsim%sFirst_Stage_RF_Data_Beamformed.mat',loadpath,filesep,filesep);
            load(str,'xmt');
        end
                
        [rcv.scanline_ref_point rcv.scanline_direction rcv.scanline_angle] = calc_scanline_position(useCaseParams);
        rcv.no_lines = size(rcv.scanline_ref_point,1);
        
        



        if(~isempty(phase_profile))
            VS = CalculateVSPosition(rcv,phase_profile,useCaseParams);
%             figure,plot(rcv.focus_point(:,1),rcv.focus_point(:,3))
%             hold on
%             plot(VS(:,1),VS(:,2))
            
            for k = 1:rcv.no_lines
                rcv.focus_point(k,[1 3]) = VS(k,:);
                
            end
        end
        
        bft3.apod.f_number        = repmat(useCaseParams.bfrcvparams(1).rcvfnum,rcv.no_lines,1) ;
        % If usecase is from ProFocus scanner update stopdepth to fit data
        % length
%         if(strcmpi('sasbprofocus',type))
%             str = sprintf('%sFirst_Stage_RF_Data_Beamformed\\%s\\First_Stage_RF_Data_Beamformed.mat',loadpath,type);            
%             load(str,'RFdata')
%             useCaseParams.scanparams(1).scanareadscr.stopdepth = (size(RFdata,1)-1)*useCaseParams.scanparams(1).c_sound/useCaseParams.bfrcvparams(1).smpfreq/2-useCaseParams.scanparams(1).scanareadscr.startdepth;
%         end

    case {'drf'}
        % remove any old directory
        if(isdir([loadpath 'First_Stage_RF_Data_Beamformed\drfsim' ]))
           rmdir([loadpath 'First_Stage_RF_Data_Beamformed\drfsim'],'s')
        end
        mkdir([loadpath 'First_Stage_RF_Data_Beamformed\drfsim'])

        rcv = Generate_Scanlines(useCaseParams, no_active_elements, symmetric, 'receive',debug); 
        rcv.no_lines = size(rcv.scanline_ref_point,1);
        bft3.apod.f_number    = useCaseParams.bfrcvparams(1).rcvfnum;
    otherwise
        error('This method is not implemented')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize bft3 and set sampling frequency and speed of sound
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bft3_system('fs',useCaseParams.bfrcvparams(1).smpfreq,'c',useCaseParams.scanparams(1).c_sound);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% bft3 - line
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bft3.line.dr              = useCaseParams.scanparams(1).c_sound/useCaseParams.bfrcvparams(1).smpfreq/2;

bft3.line.direction       = rcv.scanline_direction;

bft3.line.angle           = rcv.scanline_angle;
bft3.line.length_lines    = useCaseParams.scanparams(1).stopdepthq-...
                            useCaseParams.scanparams(1).startdepthq;

% Calculate layer thickness
aclayerthickness = useCaseParams.acmodparams(1).layerthickness1 + ...
                   useCaseParams.acmodparams(1).layerthickness2 + ...
                   useCaseParams.acmodparams(1).layerthickness3;
% Calculate origin which is relative to the receiving scan line reference
% point
% Calculate line origin which is normaly at the transducer surface and not
% on the element surface
bft3.line.origin          = rcv.scanline_ref_point + bft3.line.direction*aclayerthickness;
% In bft3 the line must start where the first sample is wanted.
% So we correct the origin with startdepthq along the direction of the line
bft3.line.origin          = bft3.line.origin +  bft3.line.direction*useCaseParams.scanparams(1).startdepthq;

bft3.line.tstart          = 2*(aclayerthickness+useCaseParams.scanparams(1).startdepthq)./useCaseParams.scanparams(1).c_sound;
bft3.line.N_points_in_line = floor(bft3.line.length_lines/bft3.line.dr);


bf_image_bft3 = zeros(bft3.line.N_points_in_line,rcv.no_lines);

% --------------------------
% bft3 Apertures (TRM, RCV)
% --------------------------
% define transducer
if(isempty(scanlines))
    scanlines = 1:rcv.no_lines;
end       

if(display == 1)
    fprintf('BFT III - Starting to beamform %d lines\n',length(scanlines));
end

switch(lower(type))
    case {'sasb first stage','srf'}
         TransducerDescription = transducer('id','8820e',...
                           'aclayerthickness',useCaseParams.acmodparams(1).layerthickness1 + ...
                                              useCaseParams.acmodparams(1).layerthickness2 + ...
                                              useCaseParams.acmodparams(1).layerthickness3,...
                           'ROC',useCaseParams.acmodparams(1).shellradius,...
                           'nr_elements_x',useCaseParams.acmodparams(1).elements,...
                           'elevation_focus',0.08); 
                       
        rcv_aperture = bft3_aperture('pos',TransducerDescription.element_positions);
        emit_aperture = bft3_aperture('pos',TransducerDescription.element_positions);
        
        for line_no = scanlines
            if(rem(line_no,ceil(length(scanlines)/10))== 0 && display == 1)
                fprintf('BFT III: SRF line no: %d of %d\n',line_no,length(scanlines));
            end
            
            % ------------------------
            % bft3 Apodization
            % ------------------------
            emit_apodizations =  bft3_apodization(emit_aperture, ...
                                    'ref',xmt.scanline_ref_point(line_no,:),  ...
                                    'distances', 0, ...
                                    'values',ones(useCaseParams.acmodparams(1).elements,1));
%             emit_apodizations.manual = true;
            emit_apodizations.fixed = true;
            
            rcv_apodizations =  bft3_apodization(rcv_aperture, ...
                                    'ref',rcv.scanline_ref_point(line_no,:),  ...
                                    'distances', 0, ...
                                    'values',rcv.apo(line_no,:)');
%             rcv_apodizations.manual = true;
             rcv_apodizations.fixed = true;
            % ------------------------
            % bft3 Line
            % ------------------------
            bft_lines = bft3_line(bft3.line.origin(line_no,:),...
                                  bft3.line.direction(line_no,:),...
                                  bft3.line.dr, ...
                                  bft3.line.length_lines);                    
                              
            % ------------------------
            % bft3 Image
            % ------------------------
            bft_image = bft3_image(emit_aperture, rcv_aperture, ...
                                   emit_apodizations, rcv_apodizations,...
                                   bft_lines);
            bft_image.interp = 'spline';
            bft_image.nthreads = int32(1);

            % ------------------------
            % bft3 Aperture focus
            % ------------------------
            emit_aperture.center_focus = xmt.scanline_ref_point(line_no,:); 
            emit_aperture.focus = xmt.focus_point(line_no,:);
            rcv_aperture.focus = rcv.focus_point(line_no,:);

            % ------------------------
            % Read channel data
            % ------------------------
            str = sprintf('%sRF_Data%sRF_Data_line_nr_%d.mat',loadpath,filesep,line_no);            
            pp = open(str);

            % ------------------------
            % bft3 Beamform
            % ------------------------
            temp = pp.RFdata;

            my_image = bft_image.beamform(temp, pp.tstart, uint32(1));
            bf_image_bft3(:,line_no) = my_image(:);
            
            figure(1)
            imagesc(bf_image_bft3)
            drawnow

        end  
      
    case {'sasbprofocus','sasbsim_new'}
        % ---------------------------------------------------------------
        % Create the emitting and receiving aperture
        % In SASB 2nd stage there is only one emitting and receiving
        % aperture for each first stage scan line. The position is at the
        % VS of the first stage scan line.
        % ---------------------------------------------------------------
        emit_aperture = bft3_aperture('pos', [0 0 0]); 
        rcv_aperture = bft3_aperture('pos', [0 0 0]);
            
        % -------------------------------------------------------------
        % Define arrays for image lines, rcv, and transmit apodizations
        % -------------------------------------------------------------
        bft_lines = [];
        emit_apodizations = [];
        rcv_apodizations = [];
        
        % -------------------------------------------------------------
        % For each image line in the final image, create an apodization 
        % The variable scanlines holds the index of the image lines to
        % beamform
        % -------------------------------------------------------------
        for line_no = scanlines
            % ---------------------------------------------------------
            % Transmit apodization
            % The transmit apodization is dynamic and ensures that only
            % image points that have been insonified by the first stage
            % scan line is included in the second stage image.
            % The variable bft3.apod holds info about f-number and
            % apodization type Hamming / Tukey / Rectwin
            % ---------------------------------------------------------
            tmp =  bft3_apodization(emit_aperture);
            tmp.fixed = false;
            tmp.dynamic = true;
            tmp.f = bft3.apod.f_number(line_no);
            tmp.window = 'Tukey';
            tmp.window = 'Hamming';
%              tmp.window = 'Rectwin';
            tmp.window_parameter  = 0.2;
            emit_apodizations = [emit_apodizations tmp];

            % ---------------------------------------------------------
            % Receive apodization
            % The receive apodization is dynamic and ensures that first 
            % stage scan lines are apodized
            % The variable bft3.apod holds info about f-number and
            % apodization type Hamming / Tukey / Rectwin
            % ---------------------------------------------------------
            tmp =  bft3_apodization(rcv_aperture);
            tmp.fixed = false;
            tmp.dynamic = true;
            tmp.f = bft3.apod.f_number(line_no);
            tmp.window = bft3.apod.type;
            rcv_apodizations = [rcv_apodizations tmp];
                
            % ------------------------
            % Image lines
            % The final image consist of a number of image lines.
            % The variable bft3.line holds information about origin,
            % direction, sample distance and line length
            % ------------------------
            tmp = bft3_line(bft3.line.origin(line_no,:),...
                            bft3.line.direction(line_no,:),...
                            bft3.line.dr, ...
                            bft3.line.length_lines);
            bft_lines = [bft_lines tmp];
        end
        
        % -------------------------------------------------
        % Open the file with first stage rf-data and tstart
        % -------------------------------------------------
        if(isempty(data))
            str = sprintf('%sFirst_Stage_RF_Data_Beamformed%s%s%sFirst_Stage_RF_Data_Beamformed.mat',loadpath,filesep,'sasbsim',filesep);            
            pp = open(str);   
            tstart = pp.tstart;
            data = pp.RFdata;
        end
        % Generate scaling matrix
        if(generate_scaling == true)
            data = ones(size(data));
        end
        
        % ------------------------
        % bft3 Beamform
        % ------------------------
%          bf_temp = zeros(size(bf_image_bft3,1),size(bf_image_bft3,2),xmt.no_lines);
        for line_no = 1:xmt.no_lines   
            % ----------------
            % Display progress
            % ----------------
            if(rem(line_no,ceil(xmt.no_lines/10))== 0 && display == 1)
                fprintf('BFT III: SASB line no: %d of %d\n',line_no,xmt.no_lines);
            end
            % -------------------------------------------------------
            % Calculate the orientation of the emitting aperture. The
            % orientation is equal to the orientation of the VS.
            % The orientation is important to update to ensure a correct 
            % transmit apodization
            % -------------------------------------------------------
            % Determine orientation vector
            orien_vec = xmt.focus_point(line_no,:)-xmt.scanline_ref_point(line_no,:);

            % Convert cartesian coordinates to spherical
            r = sqrt(sum(orien_vec.^2));
            if r == 0
                theta = 0;
            else
                theta = acos(orien_vec(3)/r); % Inclination
            end
            phi = atan2(orien_vec(2), orien_vec(1)); % Azimuth

            % Rotate so the VS is in the YZ-plane
            rot1 = phi-pi/2;

            % Rotate so the VS is on the Z-axis
            rot2 = -theta;

            % Rotational symmetry, so no need to rotate the aperture back
            rot3 = 0;

            % Rotate he emitting aperture
            emit_aperture.orientation = [rot1 rot2 rot3];
            rcv_aperture.orientation = [rot1 rot2 rot3]; % This is new and not understod
            % -------------------------------------------------------
            % Set the reference of the transmit apodizations to ensure
            % correct transmit apodization. (determines which image points
            % to include in the final image)
            % -------------------------------------------------------
            for r_idx = 1:length(scanlines)
                emit_apodizations(r_idx).ref = xmt.focus_point(line_no,:);
                rcv_apodizations(r_idx).ref = xmt.focus_point(line_no,:);
            end
               
            % ----------------------------------------------------------
            % Update the position, center focus and focus of the emitting
            % and receiving aperture
            % ----------------------------------------------------------
            emit_aperture.pos = xmt.focus_point(line_no,:);
            emit_aperture.center_focus = xmt.scanline_ref_point(line_no,:);
            emit_aperture.focus = xmt.focus_point(line_no,:);
            
            rcv_aperture.pos = xmt.scanline_ref_point(line_no,:);
            rcv_aperture.focus = xmt.focus_point(line_no,:);
            
            % ---------------
            % Setup the image
            % ---------------
            bft_image = bft3_image(emit_aperture, rcv_aperture, ...
                       emit_apodizations, rcv_apodizations,...
                       bft_lines);
                   
            bft_image.interp = 'spline';
            
            if(ispc)
                bft_image.nthreads = int32(System.Environment.ProcessorCount);
            else
                bft_image.nthreads = feature('numcores');
            end
            
            % ------------------
            % Beamform the image
            % ------------------
            if(isreal(data));
            my_image = bft_image.beamform(hilbert(data(:,line_no)), ... % Input data
                                          tstart(line_no), ...           % Start time of data [sec]
                                          uint32(1));                        % Index of emitting aperture
            else
            my_image = bft_image.beamform(data(:,line_no), ... % Input data
                                          tstart(line_no), ...           % Start time of data [sec]
                                          uint32(1));                        % Index of emitting aperture
            end
            bf_image_bft3(:,scanlines) = bf_image_bft3(:,scanlines) + my_image; 
%             bf_temp(:,:,line_no) = my_image;        
            
%             
%             figure(200)
%             imagesc(20*log10(abs(bf_image_bft3)))
%             drawnow
%      
%             figure(30)
%             imagesc(20*log10(abs(my_image)))
%             title(['line_no: ' num2str(line_no)])
%             drawnow
            
%             figure
%             plot3(xmt.focus_point(:,1),xmt.focus_point(:,2),xmt.focus_point(:,3),'.')
%             hold on, plot3(xmt.scanline_ref_point(:,1),xmt.scanline_ref_point(:,2),xmt.scanline_ref_point(:,3),'.r')
%             plot3(bft3.line.origin(:,1),bft3.line.origin(:,2),bft3.line.origin(:,3),'.k')
%             




        end
        
        
     % Generate scaling matrix
        if(generate_scaling == true)
           
            scaling = abs(bf_image_bft3);

         
            scaling = 1./scaling;
            scaling(isinf(scaling)) = 0;
            scaling(scaling > 1) = 1;
        end

    %%
    case {'sasbprofocus','sasbsim'}
        emit_aperture = bft3_aperture('pos',[0 0 0]); 
        rcv_aperture = bft3_aperture('pos',[0 0 0]);

        bft_lines = [];
        emit_apodizations = [];
        rcv_apodizations = [];
        
        % ------------------------
        % For each image line in the final image, create an apodization 
        % ------------------------
        
        for line_no = scanlines 
            % ------------------------
            % bft3 Apodization
            % ------------------------

            tmp =  bft3_apodization(emit_aperture, ...
                                    'ref',[0,0,0],  ... % [0,0,0]
                                    'distances', 0, ...
                                    'values',1);

%             tmp.manual = false;
            tmp.fixed = false;

            emit_apodizations = [emit_apodizations tmp];

            tmp =  bft3_apodization(rcv_aperture, ...
                                    'ref',[0,0,0],  ...
                                    'distances', 0, ...
                                    'values',1);
%             tmp.manual = false;
            tmp.fixed = false;
            tmp.dynamic = true;
            tmp.window = bft3.apod.type;
            tmp.f = bft3.apod.f_number(line_no);
            rcv_apodizations = [rcv_apodizations tmp];
            
            % ------------------------
            % bft3 Line
            % ------------------------
            tmp = bft3_line(bft3.line.origin(line_no,:),...
                            bft3.line.direction(line_no,:),...
                            bft3.line.dr, ...
                            bft3.line.length_lines);
            bft_lines = [bft_lines tmp];
        end
        
        bft_image = bft3_image(emit_aperture, rcv_aperture, ...
                       emit_apodizations, rcv_apodizations,...
                       bft_lines);
        bft_image.interp = 'spline';
        if(ispc)
            bft_image.nthreads = int32(System.Environment.ProcessorCount);
        else
            bft_image.nthreads = feature('numcores');
        end
        
        % ------------------------
        % Read channel data
        % ------------------------
        if(isempty(data))
            str = sprintf('%sFirst_Stage_RF_Data_Beamformed%s%s%sFirst_Stage_RF_Data_Beamformed.mat',loadpath,filesep,'sasbsim',filesep);            
            str = '/home/mah/Matlab code/SASB_test/first_stage.mat';
            pp = open(str);   
            pp.bft3.line.tstart = 2*(aclayerthickness+useCaseParams1.scanparams(1).startdepthq)./useCaseParams.scanparams(1).c_sound
            tstart = pp.bft3.line.tstart;
            data = pp.RFdata;
        end
%         str = sprintf('%sFirst_Stage_RF_Data_Beamformed%s%s%sFirst_Stage_RF_Data_Beamformed.mat',loadpath,filesep,lower(type),filesep);            
%         pp = open(str);   
%         
        % ------------------------
        % bft3 Beamform
        % ------------------------
%        vec = zeros(500,rcv.no_lines);
        for line_no = 1:xmt.no_lines    % For each
            if(rem(line_no,ceil(xmt.no_lines/10))== 0 && display == 1)
                fprintf('BFT III: SASB line no: %d of %d\n',line_no,xmt.no_lines);
            end
            % -------------------------------------------------------
            % Calculate the orientation of the emitting aperture. The
            % orientation is equal to the orientation of the VS.
            % The orientation is important to update to ensure a correct 
            % transmit apodization
            % -------------------------------------------------------
            % Determine orientation vector
            orien_vec = xmt.focus_point(line_no,:)-xmt.scanline_ref_point(line_no,:);

            % Convert cartesian coordinates to spherical
            r = sqrt(sum(orien_vec.^2));
            if r == 0
                theta = 0;
            else
                theta = acos(orien_vec(3)/r); % Inclination
            end
            phi = atan2(orien_vec(2), orien_vec(1)); % Azimuth

            % Rotate so the VS is in the YZ-plane
            rot1 = phi-pi/2;

            % Rotate so the VS is on the Z-axis
            rot2 = -theta;

            % Rotational symmetry, so no need to rotate the aperture back
            rot3 = 0;

            % Rotate he emitting aperture
            rcv_aperture.orientation = [rot1 rot2 rot3]; % This is new and not understod
            % -------------------------------------------------------
            % Set the reference of the transmit apodizations to ensure
            % correct transmit apodization. (determines which image points
            % to include in the final image)
            % -------------------------------------------------------
            for r_idx = 1:rcv.no_lines
                rcv_apodizations(r_idx).ref = xmt.focus_point(line_no,:);
            end
            
            emit_aperture.pos = [0 0 0]; % is only used if transmit apodization is used and no virtual source
            emit_aperture.center_focus = xmt.scanline_ref_point(line_no,:);
            emit_aperture.focus = xmt.focus_point(line_no,:);
            
            rcv_aperture.pos = xmt.scanline_ref_point(line_no,:);
            rcv_aperture.focus = xmt.focus_point(line_no,:);
            
            
%             tstart = pp.bft3.line.tstart(:,line_no);
            my_image = bft_image.beamform(hilbert(data(:,line_no)),...
                                          tstart,uint32(1));

            apo = FixedApod_ver2(xmt, bft3, line_no);

%            test = my_image.*apo(:,scanlines);
%            vec(:,line_no) = test(2300-249:2300+250,110);     
try
    bf_image_bft3(:,scanlines) = bf_image_bft3(:,scanlines) + my_image.*apo(:,scanlines); 

    
            figure(31)
%             imagesc(abs(hilbert(my_image)))
            subplot(1,2,1)
            imagesc(apo)
            subplot(1,2,2)
            imagesc(abs(bf_image_bft3))
            drawnow
%     imagesc(vec)
%             drawnow

        

%         figure(200)
%         imagesc(abs(bf_image_bft3))
%         
%         figure(3)
%         imagesc(abs(my_image))
%         drawnow


catch k
    keyboard
end

%     Data  = my_image.*apo(:,scanlines);
%     Data_env =  abs(hilbert(my_image.*apo(:,scanlines)));
%     save(['Data\line_no_' num2str(line_no) '_LR'],'Data')
%          video_image(line_no,:,:) = abs(hilbert(my_image.*apo(:,scanlines)));
% %             % VISUALIZE
% %             % axis(1) is the transducer setup at 1st stage line
% %             % axis(2) is the low resolution image
% %             % axis(3) is the high resolution image
%             if(exist('f1') == 0)
%                 f1 = figure(1);
%                 clf
%                 set(f1,'position',[491         592        1235         482])
%                 % Create axes
% %                 ax1 = axes('Parent',f1,'Position',[-0.0183    0.1000    0.2500    0.8000]);
%                 ax2 = axes('Parent',f1,'Position',[0.2195    0.1000    0.1400    0.8000]);
% %                 ax3 = axes('Parent',f1,'Position',[0.4102    0.1000    0.2500    0.8000]);
% %                 ax4 = axes('Parent',f1,'Position',[0.7100    0.1000    0.2500    0.8000]);
%             end
%    if(line_no > 100)         
%             Data_env =  abs(hilbert(my_image.*apo(:,scanlines)));
%             Data_env = Data_env./max(Data_env(:));
%             [Data_lg,Data_gray] = logCompress(Data_env,'Linear', 8, 60);
% 
%          
%             
%             % Scan convert
% [frame, WindowTissueX , WindowTissueY] = ScanConvert(Data_gray,...           % Data
%     useCaseParams.scanparams(1).scanareadscr.startlineorigin.x,...           % StartLineX
%     useCaseParams.scanparams(1).scanareadscr.startlineorigin.y,...           % StartLineY
%     useCaseParams.scanparams(1).scanareadscr.startlineangle,...              % StartLineAngle
%     useCaseParams.scanparams(1).scanareadscr.startdepth,...                  % StartDepth 
%     useCaseParams.scanparams(1).scanareadscr.stoplineorigin.x,...            % StopLineX
%     useCaseParams.scanparams(1).scanareadscr.stoplineorigin.y,...            % StopLineY
%     useCaseParams.scanparams(1).scanareadscr.stoplineangle,...               % StopLineAngle
%     useCaseParams.scanparams(1).scanareadscr.startdepth+bft3.line.length_lines,...% StopDepth
%     useCaseParams.scanparams(1).windowtissueq.x_tismin,...                   % WindowTissueXMin 
%     useCaseParams.scanparams(1).windowtissueq.y_tismin,...                   % WindowTissueYMin 
%     useCaseParams.scanparams(1).windowtissueq.x_tismax,...                   % WindowTissueXMax
%     useCaseParams.scanparams(1).windowtissueq.y_tismax);                     % WindowTissueYMax
%      imagesc(frame)
%      drawnow
%      figure(1)
%    end 
            
%             axes(ax1);
%             plot(TransducerDescription.element_positions(:,1)*1000,TransducerDescription.element_positions(:,3)*1000,'s','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',12)
%             hold on
%             p1 = plot(TransducerDescription.element_positions(rcv.valid(line_no,:) == 1,1)*1000,TransducerDescription.element_positions(rcv.valid(line_no,:) == 1,3)*1000,'s','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',12);
%             p2 = plot(0,80,'x','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',12);
%             
%             % calculate insonified area
%             element_nr = find(rcv.valid(line_no,:)== 1,1,'first');
%             delta_x = rcv.focus_point(line_no,1)-TransducerDescription.element_positions(element_nr,1);
%             delta_y = rcv.focus_point(line_no,3)-TransducerDescription.element_positions(element_nr,3);
%             h = 0.1;
%             A_1 = atan(delta_x/delta_y);
%             a = sin(A_1)*h;
%             b = cos(A_1)*h;
%             x_new = TransducerDescription.element_positions(element_nr,1)+a;
%             y_new = TransducerDescription.element_positions(element_nr,3)+b;
%             
%             plot([TransducerDescription.element_positions(element_nr,1) x_new]*1000,[TransducerDescription.element_positions(element_nr,3) y_new]*1000,'--','color',[0 0 0],'linewidth',2,'DisplayName','Active scanline')
%             
%             element_nr = find(rcv.valid(line_no,:)== 1,1,'last');
%             delta_x = rcv.focus_point(line_no,1)-TransducerDescription.element_positions(element_nr,1);
%             delta_y = rcv.focus_point(line_no,3)-TransducerDescription.element_positions(element_nr,3);
%             h = 0.1;
%             A_2 = atan(delta_x/delta_y);
%             a = sin(A_2)*h;
%             b = cos(A_2)*h;
%             x_new = TransducerDescription.element_positions(element_nr,1)+a;
%             y_new = TransducerDescription.element_positions(element_nr,3)+b;
%             plot([TransducerDescription.element_positions(element_nr,1) x_new]*1000,[TransducerDescription.element_positions(element_nr,3) y_new]*1000,'--','color',[0 0 0],'linewidth',2,'DisplayName','Active scanline')
%            
%             % plot half circles
%             NOP = 40;
%             center = [rcv.focus_point(line_no,1) rcv.focus_point(line_no,3)];
%             dist1 = sqrt((center(1)-TransducerDescription.element_positions(find(rcv.valid(line_no,:)== 1,1,'first'),1))^2 + ...
%                          (center(2)-TransducerDescription.element_positions(find(rcv.valid(line_no,:)== 1,1,'first'),3))^2);
%             dist2 = sqrt((center(1)-TransducerDescription.element_positions(find(rcv.valid(line_no,:)== 1,1,'last'),1))^2 + ...
%                          (center(2)-TransducerDescription.element_positions(find(rcv.valid(line_no,:)== 1,1,'last'),3))^2);
%                     
%             for k = 1:30
%                 radius = 0.003*k;
%             
%                 
%                 % draw top half circle
%                 if(radius < min(dist1,dist2))      
%                     THETA = linspace((2*pi-pi/2)-A_1,(2*pi-pi/2)-A_2,NOP);
%                     RHO=ones(1,NOP)*radius;
%                     [X,Y] = pol2cart(THETA,RHO);
%                     X=X+center(1);
%                     Y=Y+center(2);
%                     hold on
%                     plot(X*1000,Y*1000,'-b','linewidth',2);
%                 end
%                 THETA = linspace((pi/2)-A_1,(pi/2)-A_2,NOP);
%                 RHO=ones(1,NOP)*radius;
%                 [X,Y] = pol2cart(THETA,RHO);
%                 X=X+center(1);
%                 Y=Y+center(2);
%                 hold on
%                 p2 = plot(X*1000,Y*1000,'-b','linewidth',2);
%             end
%             axis ij
%             axis([rcv.scanline_ref_point(line_no,1)-32*TransducerDescription.pitch rcv.scanline_ref_point(line_no,1)+32*TransducerDescription.pitch -0.001 0.05]*1000)
%             title(['Scanline nr.: ' num2str(line_no)],'fontsize',14)
%             xlabel('Lateral Position [mm]','fontsize',12)
%             ylabel('Axial Position [mm]','fontsize',12)
%             hold off
%             axis image
% %             axis([min(rcv.focus_point(:,1)) max(rcv.focus_point(:,1)) -0.02 max(0.0985,max(rcv.focus))]*1000)
%             axis([-100 100 -20 140 ])
%             
%             axes(ax2);
%             z = (0:bft3.line.N_points_in_line-1).*bft3.line.dr+bft3.line.origin(1,3);
%             data = pp.RFdata(:,line_no);
%             data = data./max(abs(data));
%             plot(data,z*1000)
%             title('RF data','fontsize',14)
%             xlabel('Normalized amplitude','fontsize',12)
%             ylabel('Position [mm]','fontsize',12)
% %             axis image
%             axis ij
% %             axis([min(data) max(data) 75 max(0.0985,rcv.focus)*1000])
%             axis([-1 1 75 85])
%             
%             axes(ax3)
%             z = (0:bft3.line.N_points_in_line-1).*bft3.line.dr+bft3.line.origin(1,3);
%             data = my_image.*FixedApod(rcv, bft3, line_no);
%             data = 20*log10(abs(hilbert(data./max(data(:)))));
%             data(data < -60) = -60;
%             imagesc(bft3.line.origin(:,1)*1000,...
%                     z*1000,...
%                     data)
% 
%             axis ij
% %             axis([bft3.line.origin(end,1) bft3.line.origin(1,1) 0 0.085]*1000)
%             axis([-100 100 -20 140 ])
%             title('Low resolution image','fontsize',14)
%             xlabel('Lateral Position [mm]','fontsize',12)
%             ylabel('Axial Position [mm]','fontsize',12)
%             
%             axes(ax4)
%             data = 20*log10(abs(hilbert(bf_image_bft3./max(bf_image_bft3(:)))));
%             data(data < -60) = -60;           
%             imagesc(bft3.line.origin(:,1)*1000,...
%                     z*1000,...
%                     data)
%             
%             axis ij
% %             axis([bft3.line.origin(1,1) bft3.line.origin(end,1) 0.075 0.085]*1000)
%             axis([-100 100 -20 140 ])
%             title('High resolution image','fontsize',14)
%             xlabel('Lateral Position [mm]','fontsize',12)
%             ylabel('Axial Position [mm]','fontsize',12)
% 
%             drawnow
%             
%             Fra(line_no) = getframe(f1);
%             
            
        end

%%  Save frames to avi file     
%         clear bob
%         for line_no = 1:190
%            bob(line_no).cdata = Fra(line_no).cdata(20:end,1:1200,:);
%            bob(line_no).colormap = colormap;
%         end
%         
%         movie2avi(bob, 'SASB1.avi', 'compression', 'none','fps',5);
    case {'drf'}
            TransducerDescription = transducer('id','8670',...
                           'aclayerthickness',useCaseParams.acmodparams(1).layerthickness1 + ...
                                              useCaseParams.acmodparams(1).layerthickness2 + ...
                                              useCaseParams.acmodparams(1).layerthickness3,...
                           'ROC',useCaseParams.acmodparams(1).shellradius,...
                           'nr_elements_x',useCaseParams.acmodparams(1).elements,...
                           'nr_elements_y',useCaseParams.acmodparams(1).elements,...
                           'elevation_focus',0); % elevation focus set to zeros for 2d arrays on test bassis

                       
        rcv_aperture = bft3_aperture('pos',TransducerDescription.element_positions);
        emit_aperture = bft3_aperture('pos',TransducerDescription.element_positions);
        
        for line_no = scanlines
             if(rem(line_no,ceil(length(scanlines)/10))== 0 && display == 1)
                fprintf('BFT III: DRF line no: %d of %d\n',line_no,length(scanlines));
             end
            % ------------------------
            % bft3 Apodization
            % ------------------------
            emit_apodizations =  bft3_apodization(emit_aperture, ...
                                    'ref',rcv.scanline_ref_point(line_no,:),  ...
                                    'distances', 0, ...
                                    'values',ones(useCaseParams.acmodparams(1).elements^2,1));
            emit_apodizations.fixed = false;
            
            
%             rcv_apodizations = bft3_apodization(rcv_aperture, ...
%                                     'ref',rcv.scanline_ref_point(line_no,:),  ...
%                                     'distances', 0, ...
%                                     'values',reshape(rcv.valid(line_no,:),[],1));
            rcv_apodizations = bft3_apodization(rcv_aperture, ...
                                    'ref',rcv.scanline_ref_point(line_no,:),  ...
                                    'distances', 0, ...
                                    'values',ones(useCaseParams.acmodparams(1).elements^2,1));
            rcv_apodizations.fixed = true;
            rcv_apodizations.dynamic = true;
            rcv_apodizations.window = bft3.apod.type;
            rcv_apodizations.f = bft3.apod.f_number;
            % ------------------------
            % bft3 Line
            % ------------------------
            bft_lines = bft3_line(bft3.line.origin(line_no,:),...
                                  bft3.line.direction(line_no,:),...
                                  bft3.line.dr, ...
                                  bft3.line.length_lines);                    
                              
            % ------------------------
            % bft3 Image
            % ------------------------
            bft_image = bft3_image(emit_aperture, rcv_aperture, ...
                                   emit_apodizations, rcv_apodizations,...
                                   bft_lines);
            bft_image.interp = 'spline';
            bft_image.nthreads = int32(1);

            % ------------------------
            % bft3 Aperture focus
            % ------------------------
            emit_aperture.center_focus = rcv.scanline_ref_point(line_no,:); % data nulpunkt sample = t0
            emit_aperture.focus = rcv.focus_point(line_no,:);
            rcv_aperture.focus = rcv.focus_point(line_no,:);

            % ------------------------
            % Read channel data
            % ------------------------
        
            str = sprintf('%sRF_Data\\RF_Data_line_nr_%d.mat',loadpath,line_no);            
            pp = open(str);

            % ------------------------
            % bft3 Beamform
            % ------------------------
            my_image = bft_image.beamform(pp.RFdata, pp.tstart, uint32(1));
            bf_image_bft3(:,line_no) = my_image(:);
        end
%%        
end


try
switch(lower(type))
    case {'sasb first stage','srf'}
        if(~isdir([loadpath 'First_Stage_RF_Data_Beamformed' filesep 'sasbsim']))
            mkdir([loadpath 'First_Stage_RF_Data_Beamformed' filesep 'sasbsim'])
        end
        str = sprintf('%sFirst_Stage_RF_Data_Beamformed%ssasbsim%sFirst_Stage_RF_Data_Beamformed.mat',loadpath,filesep,filesep);
        RFdata = bf_image_bft3;
        save(str,'RFdata','useCaseParams','bft3','xmt','rcv');
        if(display == 1)
            fprintf('Saving first stage beamformed data\n');
        end
    case {'sasbprofocus','sasbsim_new'}        
        if(generate_scaling == true)
            save('scaling','scaling','bf_image_bft3','useCaseParams');
        else
            type = 'sasbsim';
            if(~isdir([savepath 'Second_Stage_RF_Data_Beamformed' filesep lower(type)]))
                mkdir([savepath 'Second_Stage_RF_Data_Beamformed' filesep lower(type)])
            end
            if(useCaseParams.scanparams(1).scantype == 1)
%                 keyboard
%                 % display each frame
%                 for k = 1:10:size(bf_image_bft3,1) 
%                     d = abs(bf_image_bft3(k,:));
%                     m = abs(bf_image_bft3(:));
%                     figure(1);
%                     imagesc(reshape(20*log10(d./max(m(:))),81,81))
%                     colormap gray
%                     caxis([-60 0])
%                     title(k)
%                     pause(0.1)
%                 end
%                 % display sum
%                     d = sum(abs(bf_image_bft3),1);
%                     m = max(max(sum(abs(bf_image_bft3),1)));
%                     figure(1);
%                     imagesc(reshape(20*log10(d./m),81,81))
%                     colormap gray
%                     caxis([-60 0])
%                     title('sum')
%                 RFdata = reshape(bf_image_bft3(60,:),useCaseParams.scanparams(1).stoplinenumq+1,[]);
%                 RFdata = reshape(sum(abs(bf_image_bft3),1),useCaseParams.scanparams(1).stoplinenumq+1,[]);
%                 
%                 imgxz = zeros(608,81);
%                 imgyz = zeros(608,81);
%                 for k = 1:size(bf_image_bft3,1) 
%                     d = reshape(abs(bf_image_bft3(k,:)),81,81);
%                     imgxz(k,:) = d(40,:);
%                     imgyz(k,:) = d(:,40);
%                 end
%                 m = max(imgyz(:));
%                 figure(1);
%                 imagesc(20*log10(imgyz./m))
%                 colormap gray
%                 caxis([-60 0])
%                 title('yz')
%                 figure(2);
%                 imagesc(20*log10(imgxz./m))
%                 colormap gray
%                 caxis([-60 0])
%                 title('xz')
                
                RFdata = bf_image_bft3;
                
            else
                RFdata = bf_image_bft3;
            end
            str = sprintf('%sSecond_Stage_RF_Data_Beamformed%s%s%sSecond_Stage_RF_Data_Beamformed.mat',savepath,filesep,lower(type),filesep);
            save(str,'RFdata','useCaseParams','bft3','xmt');  
            if(display == 1)
                fprintf('Saving second stage beamformed data\n');
            end
        end
        
    case {'sasbprofocus','sasbsim'}          
        if(~isdir([savepath 'Second_Stage_RF_Data_Beamformed' filesep lower(type)]))
            mkdir([savepath 'Second_Stage_RF_Data_Beamformed' filesep lower(type)])
        end
        RFdata = bf_image_bft3;
        str = sprintf('%sSecond_Stage_RF_Data_Beamformed%s%s%sSecond_Stage_RF_Data_Beamformed.mat',savepath,filesep,lower(type),filesep);
        save(str,'RFdata','useCaseParams','bft3','xmt');  
        if(display == 1)
            fprintf('Saving second stage beamformed data\n');
        end
        
    case {'drf'}
        if(~isdir([loadpath 'First_Stage_RF_Data_Beamformed\DRFsim']))
            mkdir([loadpath 'First_Stage_RF_Data_Beamformed\DRFsim'])
        end
        RFdata = bf_image_bft3;
        str = sprintf('%sFirst_Stage_RF_Data_Beamformed\\DRFsim\\First_Stage_RF_Data_Beamformed.mat',loadpath);
        save(str,'RFdata','useCaseParams');
        if(display == 1)
            fprintf('Saving dynamic receive focus beamformed data\n');
        end
    otherwise
        error('This method is not implemented')
end

catch
    keyboard      
end        
if(display == 1)
    disp('.............DONE with beamforming')
end

function VS = CalculateVSPosition(rcv,delay_profile,useCaseParams)
%%
% clear VS
for line_no = 1:192
    line_no
    x = [];
    y = [];

    delay = rcv.delays(line_no,:);
    delay = (-delay)+max(rcv.delays(line_no,:));
    
    % include delay error
    delay = delay + delay_profile;
    
    % convert from time to distance
    delay = delay*useCaseParams.scanparams(1).c_sound;
    
    for ref = 1:192
        for k = 1:192% [line_no-10 line_no line_no+10] %1:192
            if(rcv.apo(line_no,k) ~= 0)
                if(rcv.apo(line_no,ref) ~= 0)
                    d = abs(rcv.scanline_ref_point(ref,1)-rcv.scanline_ref_point(k,1));
                    R = rcv.focus(ref)+delay(ref);
                    r = rcv.focus(ref)+delay(k);

                    if(k < ref)
                        x_temp = rcv.scanline_ref_point(ref,1)+(d^2-r^2+R^2)/(2*d);
                    else
                        x_temp = rcv.scanline_ref_point(ref,1)-(d^2-r^2+R^2)/(2*d);
                    end
                    y_temp = (4*d^2*R^2-(d^2-r^2+R^2)^2)/(4*d^2);
                    y_temp = (sqrt(abs(y_temp))) - useCaseParams.acmodparams(1).layerthickness1;
                    x = [x x_temp];
                    y = [y y_temp];
                end
            end
        end

    end
    x = x(~isnan(x));
    y = y(~isnan(y));

    if(~isempty(x))
        X = [x(:), y(:)];
        if(size(X,1) > 1)
        opts = statset('Display','final');
        [cidx, ctrs] = kmeans(X, 1, 'Distance','city');
        
%         plot(X(cidx==1,1),X(cidx==1,2),'r.', ...
%              X(cidx==2,1),X(cidx==2,2),'b.', ctrs(:,1),ctrs(:,2),'kx');
        else
            ctrs = X;
        end
        VS(line_no,:) = ctrs;
    else
        VS(line_no,:) = [0 0];
    end
    
end
% keyboard
figure
hold on
set(gcf,'Position',[ 573   207   643   414])
plot(rcv.scanline_ref_point(:,1)*1000,((rcv.scanline_ref_point(:,3)+rcv.focus))*1000,'*k','linewidth',2)
plot(VS(:,1)*1000,VS(:,2)*1000,'ok','linewidth',2)
set(gca,'fontsize',14)
set(gca,'YColor',[0 0 0])
set(gca,'Units','pixels')
set(gca,'Position',[100    55   530   321])
xlabel('Lateral position [mm]','fontsize',18)
ylabel('Axial position [mm]','fontsize',18)
set(gcf,'PaperPositionMode','auto')
grid on
legend('Desired VS position','Real VS position')
drawnow