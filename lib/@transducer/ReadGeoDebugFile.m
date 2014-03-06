% Read BK Medical GeoDebugger file generated with i.e. ProFocus 2202.
%    
%   Use: This function parses a Geo Debugger file generated with Casa
%   Blanca Toolbox. The function returns a struct with all variables listed
%   in the file.
%    
%   Output:
%     Various parameters such as:
%      scan_line_table (Number X Y Angle CenterElemIdx Int_AperNum Ext_AperNum)
%      element_position_table (X Y Angle BusLine Number)
%
%   Example of use:
%     parameters = ReadGeoDebugFile(filename)
%
% Martin Hemmsen (ph.D)
% 28.03.10 (d.m.y)
% Ver. 01

function parameter = ReadGeoDebugFile(myClass_object,filename)
warning('off','MATLAB:dispatcher:InexactMatch')
warning('off','MATLAB:dispatcher:InexactCaseMatch')
parameter = [];

fid=fopen(filename);
if(fid == -1)
    parameter = [];
    fprintf('Geodebug file was not found: returning\n')
    return;
end
count = 1;
while(~feof(fid))
    tline = strtrim(fgetl(fid));

    % check for number
    if(isempty(str2num(tline(1:strfind(tline,' '))))) % not a number
        index_eq = strfind(tline,'=');
        if(~isempty(index_eq))
            index_dot = strfind(tline,'.');
            if(index_dot > index_eq)
                index_dot = [];
            end

            % retreive toptitle
            if(isempty(index_dot))
                TOPtitle = lower(tline(1:index_eq-2));
            else
                TOPtitle = lower(tline(1:index_dot-1));
            end
            % retrieve subtitle
            if(~isempty(index_dot))
                SUBtitle = lower(tline(index_dot+1:index_eq-2));
            else
                SUBtitle = [];
            end
            
            % parse value
            index_sem = strfind(tline,';');
            switch(lower(TOPtitle))
                case{'scan_line_table'}
                    fgetl(fid); % skip one line
                    
                    val = zeros(parameter.totalnumb_scanlines,7);
                    % startlinenumq starts from 0
                    for k = parameter.startlinenumq:parameter.stoplinenumq 
                        tline = strtrim(fgetl(fid));
                        val(k+1,:) = str2num(tline(1:end-1));
                    end
                case{'element_position_table'}
                    fgetl(fid); % skip one line
                    
                    val = [];
                    while(ischar(tline) && ~feof(fid))
                        tline = strtrim(fgetl(fid));
                        val = [val; str2num(tline(1:end-1))];
                    end
                otherwise
                    val = str2num(tline(index_eq+1:index_sem-1));
            end
                
            % insert field names into structure
            if(isempty(SUBtitle))
                parameter = setfield(parameter,TOPtitle,val);
            else
                parameter = setfield(parameter,TOPtitle,SUBtitle,val);
            end

        end
    end
   
     % read till end of file
    if(~ischar(tline) & tline == -1)
        disp(['Processed ' num2str(count) ' lines'])
        break
    end

    count = count +1;
end
warning('on','MATLAB:dispatcher:InexactMatch')