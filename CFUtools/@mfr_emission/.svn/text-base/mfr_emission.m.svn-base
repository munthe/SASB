
%% MFR_EMISSION
% Helper Class to set emission settings
% 
% obj = mfr_emission(key-value option pairs)
%
% Create an emission object by specifying a number of options as string-arguments pairs.
% 
% Options:
%                   c: speed of sound. --It is the only mandatory option--      (Ø)
%                type: Type of aperture. Valid options are: 
%                      vermon_32x32, dens_custom, dense_32x32, dense_64x64 and custom.
%              kerf_x: Kerf in x-dimension
%              kerf_y: Kerf in y-dimension
%             pitch_x: Pitch in x-dimension
%             pitch_y: Pitch in y-dimension
%            no_elm_x: Number of elements in x-dimension
%            no_elm_y: Number of elements in x-dimension
%              no_elm: Total number of elements
%                 pos: Position of each element. For for type='custom'. [no_elm 3]
%        no_act_elm_x: Number of active elements in x-dim             (32)
%        no_act_elm_y: Number of active elements in y-dim             (32)
%          no_act_elm: Total number active of elements
%        origin_coord: Coordinate of line origin (center of apodisation) [x y z] (mean(pos))
%          origin_elm: Element index of line origin [idx_x idx_y]  (center element)
%                      if 'origin_elm' is set, it overwrites 'origin_coord'.
%              window: Window function   ('Hanning')
%        window_alpha: Parameter to window function                 (0.75)
%          window_thr: Minimum window (APO) value   [0 1[                    (0)
%             zy_axis: Rotation around the x-axis in radians              (0rad)
%             zx_axis: Rotation around the y-axis in radians              (0rad)
%          rot_z_axis: Rotation around the z-axis in radians (has no effect) (0rad)
%             focus_r: Distance from origin to focus point [0 inf]           (Ø)
%              delays: [1024x1 double]
%                 apo: [1024x1 double]
%               apo_x: Apodization in the x-direction. cannot be used with type=custom
%               apo_y: Apodization in the y-direction. cannot be used with type=custom
%                  vs: Coordinate of the virtual source
%        elm_dead_idx: Index of the dead/very bad elements on the transducer.
%     elm_missing_idx: Index of completely missing elements (e.g., when a DAUP is missing)
%    daup_missing_idx: Index of a missing DAUP board. Sets elm_missing_idx.
% bft3_time_comp_mode: Time compensation mode sets how the BFT3 variables are
%                      calculated. It only has an influence when the VS is in fron of
%                      the array. 1 is simple, 2 is advanced and 3 is for single
%                      element emission.                                       (2)
%
%
% Read-only properties:
%             xmt_f_x: F-number calculated in the x-dimension.
%             xmt_f_y: F-number calculated in the x-dimension.
%             version: The current version number. Deprecated. Use instead:
%                      info=cfutools_get_info; info.revision
%   bft3_center_focus: Coordinate of the center_focus.
%                       -bft3_time_comp_mode=1 -> coord. of first firing elm.
%                       -bft3_time_comp_mode=2 -> coord. of center of TX aperture.
%                       -bft3_time_comp_mode=3 -> coord. of the focus point.
%          bft3_focus: Coordinate the focus point. 
%                       -bft3_time_comp_mode=1,2 -> coord. of the focus point. 
%                       -bft3_time_comp_mode=3   -> empty.
%   bft3_tstart_delta: Correction to tstart. 
%                       -bft3_time_comp_mode=1 -> 0.
%                       -bft3_time_comp_mode=2 -> Diff. in ToF. between center_focus and
%                       center of TX aperture to focus point.
%                       -bft3_time_comp_mode=3 -> ToF from center_focus to focus point.

%
% Methods:
%    plot_apo
%    plot_vs
%    plot_delay
%    plot_xdc
%             
%
% Example:  
%   ems_obj = mfr_emission('c', 1540, ...
%                          'no_act_elm',  300, ...
%                          'origin_coord', [2 2 0]/1000, ...
%                          'focus_r', 8/1000);
%  figure;
%  ems_obj.plot_vs;
% 
%
% History:
%  Version 1.0, Init version of unknown date (late 2010), M.F.Rasmussen
%  Verison 1.1, 2012-07-26, updated some months ago. Can now take an array of element positions.
%                           It now also accepts a vector representing the apodization of each element.
%  Version 1.2  2012-08-08, added type: dense_custom.
%  Version 1.3  2012-09-02, Moved apodization and element position into sub-functions
%  Version 1.4  2012-09-13, added ''elm_dead_idx'' and renamed ''missing_elm_idx'' to ''elm_missing_idx''.
%  Version 1.5  2012-09-13, bug fix.
%  Version 1.6  2012-09-14, Apodization bug fix (appeared when setting both dead and missing elements).
%  Version 1.7  2012-09-22, ''no_act_elm'' can now be used with all apertures.
%
%  Version 2.0  2012-11-20, Speed optimization. Is now ~50 times faster (apodization has been optimized).
%  Version 2.1  2013-01-06, Added a test for ''no_elm_x'' and ''no_elm_y'' when using ''custom_dense''.
%  Version 2.2  2013-03-19, Can now set missing DAUP board. Added a test for ''no_elements'' when using
%                           ''custom''. Removed ''water_temp'' as parameter.
%  Version 2.3  2013-05-XX, New convention for the dimension directions are
%                           implemented. The transducer is flipped 90 degrees. This
%                           convention has also been implemented on SARUS.
%  Version 2.4  2013-06-20, Added a 'field_2_compat_mode' -when set, the new convention
%                           from version 2.3 is rolled back. This is for easier
%                           interaction with Field 2.
%  Version 2.5  2013-08-22, Added normalisation of apodisation -ensures maximum is
%                           always 1. Added the member "first_elm" which contains the
%                           coordinate of the first element that fires. Updated
%                           plot_apo to give a prettier output. Renamed no_elements_ to
%                           no_elm_. Updated the help text. 
%  Version 2.6 2013-08-29   Speeded up by using parse_input_parameters instead of
%                           bft3_va_ar. Added the property CFUtools_info.
%  Version 2.7 2013-10-29   Moved calculation of variables used for beamforming with
%                           BFT3 from beamforming script to this function. New properties are:
%                           bft3_time_comp_mode, bft3_center_focus, bft3_focus, bft3_tstart_delta.



%

classdef mfr_emission < handle

properties (SetAccess = 'private', GetAccess = 'public')
    version       = '2.7';
    CFUtools_revision = [];
    CFUtools_info = [];
    type          = 'vermon_32x32';
    kerf_x        = 22e-6;
    kerf_y        = 22e-6;
    pitch_x       = 300e-6;
    pitch_y       = 300e-6;
    no_elm_x      = [];
    no_elm_y      = [];
    no_elm        = [];
    pos           = [];
    c;
    no_act_elm_x  = [];
    no_act_elm_y  = [];
    no_act_elm    = [];
    origin_coord  = [];   % start of line (center of apo)
    origin_elm    = [];   % start of line (center of apo)
    window        = 'Hamming';
    window_alpha  = 0.75;
    window_thr    = 0;
    zy_axis       = 0; % rotation around x-axis
    zx_axis       = 0; % rotation around y-axis
    rot_z_axis    = 0; % rotation around z-axis (no effect)(in the future might affect APO)
    focus_r       = 0;
    xmt_f_x
    xmt_f_y
    delays      % delays

    apo           = []; % apodisation
    apo_x         = [];
    apo_y         = [];
    vs          % virtual source xyz coordinate
    elm_dead_idx     = [];
    elm_missing_idx  = [];
    daup_missing_idx = [];
    first_elm        = []; %coordinate of first position that fires
    field_2_compat_mode = false;
    bft3_time_comp_mode = 2;
    bft3_center_focus   = [];
    bft3_tstart_delta   = 0;
    bft3_focus          = [];
end


properties (SetAccess = 'private', GetAccess = 'private')
    pos_x;
    pos_y;
    pos_idx_x;
    pos_idx_y;
    pos_dead;
    no_dead_rows;
    dead_row_idx;
    apo_x_idx;
    apo_y_idx;
    apo_external_set = 0;
    pos_external_set = 0;
end


methods
        function obj = mfr_emission(varargin)
        %
        % obj = mfr_emission(varargin)
        %
        % Vars: 
        %  'window' (those supported by mfr_window)                      ('Hamming')
        %  'window_alpha': parameter to window function [0 1]                 (0.75)
        %  'window_thr'  : minimum window (APO) value   [0 1[                    (0)
        %  'no_act_elm_x' : number of active elm in x-dim [1 32]           (32)
        %  'no_act_elm_y' : number of active elements in y-dim [1 32]           (32) 
        %  'origin'  : coordinate of line origin (center of apodisation)   ([0 0 0])
        %  'focus_r' : distance from origin to focus point [0 inf]              (0m)
        %  'zy_axis'  : rotation around the x-axis in radians                 (0rad)
        %  'zx_axis'  : rotation around the y-axis in radians                 (0rad)
        %  'c'           : speed of sound.                                       (Ø)
        tic;

        %% Set the standard values
        st.type          = 'vermon_32x32';
        st.c             = obj.c;   %[m/s]
        st.focus_r       = obj.focus_r; %[m]
        st.no_act_elm_x  = obj.no_act_elm_x;
        st.no_act_elm_y  = obj.no_act_elm_y;
        st.no_act_elm    = obj.no_act_elm;
        st.pitch_x       = obj.pitch_x;
        st.pitch_y       = obj.pitch_y;
        st.kerf_x        = obj.kerf_x;
        st.kerf_y        = obj.kerf_y;
        st.focus_r       = obj.focus_r;      % [m]
        st.origin_coord  = obj.origin_coord; % line start / center of active aperture
        st.origin_elm    = obj.origin_elm;   % line start / center of active aperture
        st.window        = obj.window;
        st.window_alpha  = obj.window_alpha;
        st.window_thr    = obj.window_thr;
        st.zy_axis       = obj.zy_axis;      % rotation around x-axis
        st.zx_axis       = obj.zx_axis;      % rotation around y-axis
        st.rot_z_axis    = obj.rot_z_axis;   % rotation around z-axis (only affects APO)
        st.daup_missing_idx = obj.daup_missing_idx;
        st.elm_missing_idx  = obj.elm_missing_idx;
        st.elm_dead_idx     = obj.elm_dead_idx;
        st.apo           = obj.apo;
        st.apo_x         = obj.apo_x;
        st.apo_y         = obj.apo_y;
        st.pos           = obj.pos;
        st.no_elm_x      = obj.no_elm_x;
        st.no_elm_y      = obj.no_elm_y;
        st.bft3_time_comp_mode = obj.bft3_time_comp_mode;
        st.field_2_compat_mode = obj.field_2_compat_mode;
        
        
        % overwrite standard values with values given by the user
        st = cfu_parse_input_parameters(st,varargin);

        % Copy the struct values back to the object        
        obj.type          = st.type;
        obj.c             = st.c;
        obj.focus_r       = st.focus_r;      %[m]
        obj.no_act_elm_x  = st.no_act_elm_x;
        obj.no_act_elm_y  = st.no_act_elm_y;
        obj.no_act_elm    = st.no_act_elm;
        obj.pitch_x       = st.pitch_x;
        obj.pitch_y       = st.pitch_y;
        obj.kerf_x        = st.kerf_x;
        obj.kerf_y        = st.kerf_y;
        obj.origin_coord  = st.origin_coord; % line start / center of active aperture
        obj.origin_elm    = st.origin_elm;   % line start / center of active aperture
        obj.window        = st.window;
        obj.window_alpha  = st.window_alpha;
        obj.window_thr    = st.window_thr;
        obj.zy_axis       = st.zy_axis;    % rotation around x-axis
        obj.zx_axis       = st.zx_axis;    % rotation around y-axis
        obj.rot_z_axis    = st.rot_z_axis; % rotation around z-axis (only affects APO)
        obj.daup_missing_idx = st.daup_missing_idx;
        obj.elm_missing_idx  = st.elm_missing_idx;
        obj.elm_dead_idx     = st.elm_dead_idx;
        obj.apo           = st.apo;
        obj.apo_x         = st.apo_x;
        obj.apo_y         = st.apo_y;
        obj.pos           = st.pos;
        obj.no_elm_x      = st.no_elm_x;
        obj.no_elm_y      = st.no_elm_y;
        obj.bft3_time_comp_mode = st.bft3_time_comp_mode;
        obj.field_2_compat_mode = st.field_2_compat_mode;
        
        % Set the CFUtools version -- disabled. Took too long time
        %obj.CFUtools_info     = cfu_tools_get_info;
        %obj.CFUtools_revision = obj.CFUtools_info.revision;
        
        %% --------------------------------------------------
        % Call Test function if no input given
        % ---------------------------------------------------
        if nargin == 0
            obj = obj.test();
            return;
        end

        
        
        %% --------------------------------------------------
        % Set Flags
        % --------------------------------------------------
        if ~isempty(obj.pos) 
            obj.pos_external_set = 1;
        end
        if ~isempty(obj.apo) 
            obj.apo_external_set = 1;
        end

        
        %% --------------------------------------------------
        % Input Variables Verification
        % ---------------------------------------------------
        % Sound speed test
        if isempty(obj.c)
            error(' You must set the speed of sound.\n')
        end
        if obj.window_thr < 0 || obj.window_thr >= 1,
            error('window_thr must be larger than or equal to 0, and smaller than 1.');
        end
        if ~isempty(st.no_act_elm) && (~isempty(st.no_act_elm_x) || ~isempty(st.no_act_elm_y))
            error('When ''no_act_elm'' is set, ''no_act_elm_x'' and ''no_act_elm_y'' cannot be set.');
        end
        
        
        
        
        
        %% --------------------------------------------------
        % Transducer Specific Setup (element positions)
        % --------------------------------------------------
        obj = obj.aperture_setup();
        obj = obj.element_positions();
        
        

        %% -----------------------------
        % Dead elements
        %-----------------------------
        if ~strcmp(obj.type, 'custom')
            obj = obj.remove_dead_elements();
        end

        
        %% -----------------------------
        % Apodisation 
        %-----------------------------
        if obj.apo_external_set == 0
            obj = obj.apodization();
        end
        
        
        
        % make sure the apodization length is right
        if length(obj.apo) ~= size(obj.pos,1)
            error('Length of apodization must be %i, it was: %i.',...
                  size(obj.pos,1), length(obj.apo));
        end

        
        %% -----------------------------
        % F-number
        %-----------------------------
        if (~strcmp(obj.type, 'custom'))
            if ~isempty(obj.apo_x) && ~isempty(obj.apo_y)
                active_width_x = sum(obj.apo_x > 0)/obj.pitch_x;
                active_width_y = sum(obj.apo_y > 0)/obj.pitch_y;
            elseif isempty(obj.apo_x) && isempty(obj.apo_y) && ~isempty(obj.apo)
                active_width_x = sum(obj.apo(1:32:end) > 0)/obj.pitch_x;
                active_width_y = sum(obj.apo(1:32) > 0)/obj.pitch_y;
            else
                active_width_x = abs(obj.pos_x(obj.apo_x_idx(end))-obj.pos_x(obj.apo_x_idx(1)))*abs(cos(obj.zx_axis)); 
                active_width_y = abs(obj.pos_y(obj.apo_y_idx(end))-obj.pos_y(obj.apo_y_idx(1)))*abs(cos(obj.zy_axis));
            end
            obj.xmt_f_x = obj.focus_r/active_width_x;
            obj.xmt_f_y = obj.focus_r/active_width_y;
        end

        
        %% -----------------------------
        % Focus Point (Virtual Source)
        %-----------------------------
        % TODO: Rotate by euler angles
        % Relative focus point
        obj.vs = zeros(1,3);
        obj.vs(3) = sqrt(obj.focus_r^2 / (1+tan(obj.zx_axis)^2+tan(obj.zy_axis)^2));% z[m]
        obj.vs(1) = obj.vs(3)*tan(obj.zx_axis); % x[m]
        obj.vs(2) = obj.vs(3)*tan(obj.zy_axis); % y[m]
        
        if obj.focus_r>=0  % VS translated by origin
            obj.vs = obj.vs+obj.origin_coord;
        else %if negative direction 
            obj.vs=-obj.vs+obj.origin_coord;
        end

        %% -----------------------------
        % Delays
        %-----------------------------
        obj = calc_delays(obj);
        
        
        %% -----------------------------
        % BFT-3 Variables
        %-----------------------------
        obj = calc_bft3_variables(obj);
        
        end %function
end %method
end %class def


