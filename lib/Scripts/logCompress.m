function [y_out, y_dig, y, x, xm, x0, y0, xf, yf, mu_c, alpha] = logCompress(inputData, methodSel, nOut, DR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%	inputData   linear envelope data
%   methodSel   Compression method:
%                   BK_muLaw_RTSC_VP_DRC_SEL_0
%                   BK_modif_RTSC_VP_DRC_SEL_1
%                   BK_muLaw_RTSC_VP_DRC_SEL_2
%                   BK_pure
%                   Linear
%   nOut        Number of bits for compressed output.
%   DR          Dynamic range of input data
% OUTPUT
%   y_out       Log-compressed output, dB. [-y0dB; 0]. e.g. [-48.1; 0]
%   y_dig       Log-compressed output, digital value. e.g. [0;255]
%   y           Log-compressed output, normalized to 1 e.g. [0.0039;1]
%   x           Input data after clipping
%   xm          Input resulting in max output
%   x0          Min input signal.
%   y0          Min output signal.
%   xf          Cross point for different DR (input value)
%   yf          Cross point for different DR (output value)
%   mu_c        Calculated (complex) mu (zero-cross. or newton). Only R{mu} is used in alg.
%   alpha       Modified Log parameter
% DESCRIPTION
%   -   Logarithmic compression of data with a certain input range (DR) to an output
%       with a fixed number of output bits (nOut).
%   -   Compression can be done using method of choice (method).
%   -   Output are converted to dB scale.
%
% VERSION
%   v2  2010-05-26
%       mu-law: the stated mistake in the simple approach was that xm=1 was assumed, but
%       actually the assumption is that x0/xm = 1/R. That is there is no need for the
%       expanded approach.
%   v1  2010-04-15
% AUTHOR    Jacob kortbek
%
% ################################
% SOURCE CODE:
% ################################
%
% \Engine\SourceCode\RTSC\VPCtrl\Vpcontra.c:
% --------------------------------------------
% int VP_ContrastCalcLineParams(RTSC_MidLevelMode MidLevelMode)
% {
% #if VECP_VERSION >= 0x0006
% 	int DRC_Select;
%   ... 
% 	DRC_Select = VP_Data.MidLevelModeDataTbl[MidLevelMode].DRC_Select;
%   ...
% #endif
% }
%
% \Engine\SourceCode\Common\INCLUDE\rtscdata.h
% --------------------------------------------
% RTSC_VP_Params
% Holds parameters specifying how Vector Processing is to be made.
% Description:    
%     Holds parameters specifying how Vector Processing is to be made. For MidLevelMode DUM1 and DBG1 all fields are unused.
% */
% typedef struct RTSC_VP_Params_struct 
% {
%     ...
%     RTSC_VP_DRC_Sel DRC_Select; /* Selection of the type of dynamic range 
%       compression algorithm to use. Controls B/M gray scale image contrast. 
%       Only applicable for B/M modes - set to 0 (zero) for other modes. */
%     ...
% } RTSC_VP_Params;
%
% \Engine\SourceCode\Common\INCLUDE\rtscdata.h
% --------------------------------------------
%             RTSC_VP_DRC_Sel
%             Description:
%                 Dynamic range compression algorithm selection.
%             */
%             typedef enum RTSC_VP_DRC_Sel_enum
%             {
%               RTSC_VP_DRC_SEL_0, /* My-law compression algorithm 1 (default). */
%               RTSC_VP_DRC_SEL_1, /* Modified logarithmic compression 1. */
%               RTSC_VP_DRC_SEL_2, /* My-law compression algorithm 2. */
%             } RTSC_VP_DRC_Sel;
%
% \Engine\SourceCode\RTSC\VPCtrl\Vpcontra
% --------------------------------------------
%              "contrast algorithm "
%
% \Engine\SourceCode\RTSC\VPCtrl\Vpbctrl
% --------------------------------------------
%             /* Calculate and download the contrast LUT used for dynamic range compression */
%             /* Set the contrast curve control data to the default value: */
%             m_ContrastDat[0].FixPoint = 0.5f;
%             m_ContrastDat[0].Gain = 0.0f;
%             m_ContrastDat[0].Alfa = 0.0f; /* N.A. */
%             m_ContrastDat[0].AlgSel = 0;  /* My-law */
%             m_ContrastDat[1].FixPoint = 0.25f;
%             m_ContrastDat[1].Gain = 0.0f;
%             m_ContrastDat[1].Alfa = 0.00005f;
%             m_ContrastDat[1].AlgSel = 2;  /* Modif. log. 1 */
%             m_ContrastDat[2].FixPoint = 0.35f;
%             m_ContrastDat[2].Gain = 0.0f;
%             m_ContrastDat[2].Alfa = 0.0f; /* N.A. */
%             m_ContrastDat[2].AlgSel = 0;  /* My-law */
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -----------------------------------
% Set parameters for Log-Compression
% -----------------------------------

xm                  = 1; % max input
ym                  = 1; % max output
DR_max              = 90;

if strcmp(methodSel,'BK_muLaw_RTSC_VP_DRC_SEL_0')
    yf = 0.5; % 0.5
    calcMethod = 'BK_muLaw';
elseif strcmp(methodSel,'BK_modif_RTSC_VP_DRC_SEL_1')
    yf = 0.25; % 0.25
    calcMethod = 'BK_modif';    
elseif strcmp(methodSel,'BK_muLaw_RTSC_VP_DRC_SEL_2')
    yf = 0.35; % 0.35
    calcMethod = 'BK_muLaw';    
elseif strcmp(methodSel,'BK_pure') 
    yf = 0.25; % 0.25
    calcMethod = 'BK_pure';    
elseif strcmp(methodSel,'Linear') 
    yf = 1;
    calcMethod = 'Linear';    
end


% "Modified Log"
% ----------------
alpha               = 0.00005; 
mod_expanded        = 0; % original RTSC-version or expanded (newton) method
% newton
mod_newtStart       = 10^(-30/20); % start value for xm in newton search
mod_newtN           = 100; % number of iterations


% "my-law"
% ----------------
mu_c                = []; % initialize
mu_vector           = 0:0.1:1000; % test vector
% newton
mu_newtEnable       = 1; % either newton or zero-search
mu_newtStart        = 500; % start value for mu in newton search
mu_newtN            = 6; % number of iterations
% plot
plot_uLaw           = 0;


% Calculate parameters
% ----------------
yOutMax     = (2^nOut-1); % max output level
y0          = 1/yOutMax; % min output value (BK only)
y0dB        = 20*log10(y0);
x0          = 10.0^(-DR/20.0);
R           = 10.0^(DR/20.0);
R_max       = 10.0^(DR_max/20.0);


% rename
x = inputData;


% -----------------------------------
% Log-Compression
% -----------------------------------


switch calcMethod

    case 'BK_muLaw'

        % non-linear compression with common cross point
        
        % find mu for DR=-90, xm=1 (where function-output y0 crosses zero)
        % y0 = ym*log(1+mu*abs(x0)/xm)./log(1+mu)
        x0_90dB = 10^(-DR_max/20.0);
        muFun = ym*log(1+mu_vector*x0_90dB./xm)./log(1+mu_vector)-y0;
        [dummy, idx_90dB] = find(diff(sign(muFun)) > 0);
        mu_zero_90dB = mu_vector(idx_90dB);

        % Find mu using newton raphson method
        mu_newt_90dB = mu_newtStart;
        for k = 1:mu_newtN
            mu_newt_90dB = mu_newt_90dB - ((1+mu_newt_90dB)^y0-mu_newt_90dB*x0_90dB/xm-1) / (y0*(1+mu_newt_90dB)^(y0-1)-x0_90dB/xm);
        end

        % select mu value
        if mu_newtEnable
            mu_90dB = mu_newt_90dB;
        else
            mu_90dB = mu_zero_90dB;
        end

        % Common cross-point:
        % Notice xm varies with DR when common cross-point is used.
        % Select a cross-point (xf,yf) and xm can be expressed in terms mu
        xf          = xm/mu_90dB*((1+mu_90dB).^(yf/ym)-1);
        xfdB        = 20*log10(xf);

        % find mu given R
        muFun = ym*log(1+mu_vector*x0./xm)./log(1+mu_vector)-y0;
        [dummy, idx] = find(diff(sign(muFun)) > 0);
        mu_zero = mu_vector(idx);

        % Find mu using newton raphson method. (xm=1)
        mu_newt = mu_newtStart;
        for k = 1:mu_newtN
            mu_newt = mu_newt - ((1+mu_newt)^(y0/ym)-mu_newt/R-1) ./ (y0/ym*(1+mu_newt)^(y0/ym-1)-1/R);
        end

        if plot_uLaw
            figure,hold on,grid on
            h1 = plot(mu_vector,muFun,'b.-');
            h2 = plot(mu_vector(idx),muFun(idx),'bo','markersize',10);
            if isreal(mu_newt)
                h3 = plot([mu_newt,mu_newt],[min(muFun),max(muFun)],'b--');
            end
            legend([h1(1),h2(1)],'muFun','zero-cross',2)
            if isreal(mu_newt)
                legend([h1(1),h2(1),h3(1)],'muFun','zero-cross','newton',2)                
            end
            axis tight
            xlim([-100 mu_vector(end)])
            title(['mu function and estimates. DR=',num2str(DR)])
        end

        if mu_newtEnable
            mu = mu_newt;
        else
            mu = mu_zero;
        end
        
        mu_c = mu; % calculated (complex) mu value returned from function

        if isempty(mu) || ~isreal(mu)
            displ_err(['fun_logCompress: invalid mu value. DR = ',num2str(DR)])
            parm = {DR, mu_zero,mu_newt};
            txt = {'DR','mu_zero','mu_newt'};
            displ(parm,txt,'%%')
        end

        % use only real part in algorithm
        mu = real(mu);
        
        % calculate new xm
        xm = xf.*mu./((1+mu).^(yf./ym)-1);

        % calculate new x0        
        x0 = 0; % this is done in source
        x0 = xm/R; % NOT in Source (x0 is not used, only for the output)
        
        % Clip input values to xm
        x(x > xm) = xm;

        % Do NOT Clip input values to x0 (not done i source either) 
        % (BELOW LINE MUST BE COMMENTED)
        %x(x < x0) = x0; % FOR TEST ONLY

        % calc y(x), given mu, DR, xm
        y = ym* log(1+mu.*x./xm) ./ log(1+mu);

    case 'BK_modif'
        
        % non-linear compression with common cross point

        x_min   = xm / R_max;
        xf      = (alpha + x_min) * exp(yf * log((alpha + xm)/(alpha + x_min))) - alpha;
        x0      = xm / R;
        if mod_expanded == 0
            % Original version
            % --> NO true common cross-point in yf
            % --> x0 is changed
            % --> xm/x0 = DR
            x_yf    = (alpha + x0) * exp(yf * log((alpha + xm)/(alpha + x0))) - alpha;
            r       = xf / x_yf;
            x0      = x0 * r;
            xm      = xm * r;
        elseif mod_expanded
            % JBK expanded version
            % --> True common cross-point in yf
            % --> x0 is changed
            % --> xm/x0 = DR
            
            % Find xm using newton raphson method.
            xm = mod_newtStart;
            for k = 1:mod_newtN
                xm = xm - ((((alpha+xm)/(alpha+xm/R))-((alpha+xf)/(alpha+xm/R)).^(1/yf))*(alpha+xm/R).^2)  ./ ( (alpha+xm/R)-(alpha+xm)/R+1/yf*((alpha+xf)/(alpha+xm/R)).^(1/yf-1)*(alpha+xf)/R);
            end
            x0 = xm/R;
        end

        % Clip input values to x0, xm
        x(x < x0) = x0;
        x(x > xm) = xm;

        % Calc y(x) DR Compression.
        y = log((alpha+x)/(alpha+x0)) /...
            log((alpha+xm)/(alpha+x0));

    case 'BK_pure'

        % linear compression with common cross point
        
        x_min   = xm / R_max;
        xf      = x_min * exp(yf * log(R_max)); % = x_min*R_max^yf
        x0      = xf / exp(yf * log(R)); % = xf / R^yf
        xm      = x0 * R;

        % Clip input values to x0, xm
        x(x < x0) = x0;
        x(x > xm) = xm;

        % Calc y(x) compressed output
        y = log(x/x0) / log(xm/x0);
        
    case 'Linear'

        % Linear mapping, common cross-point in (0,0) dB
        % NOT limited by nOut        

        xf      = 1; 
        yf      = 1;         
%         y0dB    = -DR;
%         y0      = 10^(y0dB/20); % min output value (BK only)   
        

        % Clip input values to x0, xm
        % x0 is the minimum value i.e. 0dB but not log compressed
        % xm is the maximum value i.e. -DR but not log compressed
        x(x < x0) = x0;
        x(x > xm) = xm;
        % Calc y(x) compressed output - y is in the range -DR to 0 dB
        y = 20*log10(x);
        % convert to be inbetween 0 and 1
        y = (y+DR)/DR;        

    otherwise
        y = [];
        disp('wrong DR compression method')
end

% Clip output values to a maximum of ym.
y(y > ym) = ym;

% Convert to n-bit digital output value
y_dig = floor(y/y0 + 0.5); % e.g. [0;255]

% Convert to "log-scale"
% y_out = y_dig*y0*(-y0dB) + y0dB ; % e.g. [-48.1;0] dB
y_out = y*DR-DR; % e.g. [-48.1;0] dB


max(y_dig(:));
min(y_dig(:));

max(y_out(:));
min(y_out(:));



