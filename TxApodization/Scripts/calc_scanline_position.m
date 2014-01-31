% Calculate position and direction of every scan line.
function [pos, dir, angle] = calc_scanline_position(useCaseParams)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%       useCaseParams
% OUTPUT
%   pos
%   dir
%   angle
% DESCRIPTION
%   Calculate position and direction of every scan line.
%   - First scan line at positive x-coordinates.
%   - Last scan line at negative x-coordinates.
%   - scan line direction is symmetric around pi/2
%   - scan line direction increases from first to last element.
% VERSION		
%   v1  2010-05-26
% AUTHOR    Martin Christian Hemmsen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

StartAngle = useCaseParams.scanparams(1).scanareadscr.startlineangle;
StopAngle = useCaseParams.scanparams(1).scanareadscr.stoplineangle;

Startlineorigin_X = useCaseParams.scanparams(1).scanareadscr.startlineorigin.x;
Stoplineorigin_X = useCaseParams.scanparams(1).scanareadscr.stoplineorigin.x;
Startlineorigin_Y = useCaseParams.scanparams(1).scanareadscr.startlineorigin.y;
Stoplineorigin_Y = useCaseParams.scanparams(1).scanareadscr.stoplineorigin.y;
StartLineNumQ = useCaseParams.scanparams(1).startlinenumq;
StopLineNumQ = useCaseParams.scanparams(1).stoplinenumq;
Aclayerthickness = useCaseParams.acmodparams(1).layerthickness1 + ...
                   useCaseParams.acmodparams(1).layerthickness2 + ...
                   useCaseParams.acmodparams(1).layerthickness3;
AROC = useCaseParams.acmodparams(1).shellradius - Aclayerthickness;
ROC = useCaseParams.acmodparams(1).shellradius;
Steeringangle = useCaseParams.scanparams(1).steeringangle;

Nlines = StopLineNumQ-StartLineNumQ+1;

switch(ROC)
    case 0
        
       switch(useCaseParams.scanparams(1).scantype)
           case(0) % 1d imaging
            pos = [linspace(Startlineorigin_X,Stoplineorigin_X,Nlines)',...
                   linspace(Startlineorigin_Y,Stoplineorigin_Y,Nlines)', ...
                   repmat(-Aclayerthickness,Nlines,1)];

            Phasedangle = linspace(-useCaseParams.scanparams(1).phasedangle,useCaseParams.scanparams(1).phasedangle,Nlines)';   
               
            angle = repmat(pi/2+Steeringangle,Nlines,1)+Phasedangle;

            % calculate unit vector of direction
%             dir = [cos(angle) zeros(Nlines,1) sin(angle)];
             % Determine scanline direction
                pasedangle = useCaseParams.scanparams(1).phasedangle;

            pasedangles = linspace(-pasedangle,pasedangle,Nlines);

            for n = 1:Nlines
                dir(n,:) = beam2cart_coord(1,0, pasedangles(n));
            end
                
           case(1) % we assume 2d imaging
                % Determine scanline ref position
                x = linspace(Startlineorigin_X,Stoplineorigin_X,Nlines);
                y = linspace(Startlineorigin_Y,Stoplineorigin_Y,Nlines);

                x = reshape(repmat(x,1,Nlines),[],1);
                y = reshape(repmat(y,Nlines,1),[],1);

                z = ones(Nlines*Nlines,1)*-Aclayerthickness;

                pos = [x y z];

                % Determine scanline angle  
                angle = repmat(pi/2+Steeringangle,Nlines^2,1);                   

                % Determine scanline direction
                pasedangle = useCaseParams.scanparams(1).phasedangle;

                pasedangles = linspace(-pasedangle,pasedangle,Nlines);

                for k = 1:Nlines
                    for n = 1:Nlines
                            dir((k-1)*Nlines+n,:) = beam2cart_coord(1, pasedangles(k), pasedangles(n));
                    end
                end
            case(2) % x-z plane
                % Determine scanline ref position
                x = linspace(Startlineorigin_X,Stoplineorigin_X,Nlines);
                y = linspace(0,0,Nlines);
                
                z = ones(1,Nlines)*-Aclayerthickness;

                pos = [x; y; z]';

                % Determine scanline angle  
                angle = repmat(pi/2+Steeringangle,Nlines,1);    
                
                % Determine scanline direction
                pasedangle = useCaseParams.scanparams(1).phasedangle;

                pasedangles = linspace(-pasedangle,pasedangle,Nlines);

                for n = 1:Nlines
                    dir(n,:) = beam2cart_coord(1, 0, pasedangles(n));
                end
            
            case(3) % y-z plane
                % Determine scanline ref position
                x = linspace(0,0,Nlines);
                y = linspace(Startlineorigin_Y,Stoplineorigin_Y,Nlines);

                z = ones(1,Nlines)*-Aclayerthickness;

                pos = [x; y; z]';

                % Determine scanline angle  
                angle = repmat(pi/2+Steeringangle,Nlines,1);    
                
                % Determine scanline direction
                pasedangle = useCaseParams.scanparams(1).phasedangle;

                pasedangles = linspace(-pasedangle,pasedangle,Nlines);

                for n = 1:Nlines
                    dir(n,:) = beam2cart_coord(1, pasedangles(n),0);
                end
            
            case(4) % we assume 2d imaging plane
                % Determine scanline ref position
                x = linspace(Startlineorigin_X,Stoplineorigin_X,Nlines);
                y = linspace(Startlineorigin_Y,Stoplineorigin_Y,Nlines);

                x = reshape(repmat(x,1,Nlines),[],1);
                y = reshape(repmat(y,Nlines,1),[],1);


                z = ones(Nlines*Nlines,1)*-Aclayerthickness;


                pos = [x y z];

                % Determine scanline angle  
                angle = repmat(pi/2+Steeringangle,Nlines,1);    

                % Determine scanline direction
                pasedangle = 0;

                for k = 1:Nlines
                    for n = 1:Nlines
                            dir((k-1)*Nlines+n,:) = [0 0 1];
                    end
                end
        end
    otherwise        
        ElementRadius       = AROC;
        ShellRadius         = ROC;
        AngleStep           = (StopAngle-StartAngle) / (Nlines-1);
        
        
        pos = zeros(Nlines,3);
        angle = zeros(Nlines,1);
        dir = zeros(Nlines,3); 
        for k = 1:Nlines
            arcx          = StartAngle+(k-1)*AngleStep;
            fRad          = ElementRadius;
            pos(k,1)      = fRad * cos(arcx);
            pos(k,3)      = fRad * sin(arcx) - ShellRadius;
%             angle(k)      = arcx+Steeringangle+Phasedangle(k);
            angle(k)      = arcx+Steeringangle;
        end 
  
        
        % calculate unit vector of direction
        dir(:,1) = cos(angle);
        dir(:,3) = sin(angle);
        

            
end







