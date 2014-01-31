%Read_Usecase - Read a usecase file (.dat) generated with the ProFocus
%               Research Interface
%
% Syntax:  [parameter] = Read_Usecase(filename)
%
% Inputs:
%    filename - filename of usecase
%
% Outputs:
%    parameter - struct with data
%
% Example: 
%   [parameter] = Read_Usecase('C:test.dat')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% Author: Martin Christian Hemmsen
% Danish Technical University and BK Medical
% email: mah@elektro.dtu.dk
% Website: http://server.elektro.dtu.dk/personal/mah/
% April 2010; Last revision: 27-April-2010
function [parameter] = Read_Usecase(USECASE_FILENAME)
warning('off','MATLAB:dispatcher:InexactMatch')
warning('off','MATLAB:dispatcher:InexactCaseMatch')
parameter = [];

fid=fopen(USECASE_FILENAME);
if(fid == -1)
    parameter = [];
    fprintf('file was not found: returning\n')
    return;
end
BF = NaN;
count = 1;
HEADLINE = [];
while 1
    tline = fgets(fid);
    
    % retreive parameter Headlines
    if(strfind(tline,'['))
        BF = str2num(tline(strfind(tline,'.')+1:strfind(tline,']')-1))+1;
        if(isempty(BF))
            HEADLINE = lower(tline(strfind(tline,'[')+1:strfind(tline,']')-1));
        else
            HEADLINE = lower(tline(strfind(tline,'[')+1:strfind(tline,'.')-1));
        end
    end
    
    % insert parameter values into parameter
    if(strfind(tline,'='))
       label = lower(tline(1:strfind(tline,'=')-1));
       if(isempty(BF))
           % special case for MacroRoutine
           if(strcmp(HEADLINE,'macroroutine'))
               index = [strfind(tline,'=') strfind(tline,',') length(tline)];
                for k = 1:4
                    parameter = setfield(parameter,HEADLINE, ['F' label],...
                    str2num(tline(index(k)+1:index(k+1)-1)));
                end
           else
                parameter = setfield(parameter,HEADLINE,label,...
                str2num(tline(strfind(tline,'=')+1:end)));
           end
       else
            index = [0 strfind(label,'.') length(label)+1];
            clear label
            for k = 1:length(index)-1
                if(isempty(str2num(tline(index(k)+1:index(k+1)-1))))
                    label{k} = lower(tline(index(k)+1:index(k+1)-1));
                else
                    label{k} = ['F' tline(index(k)+1:index(k+1)-1)];
                end
                    
            end
            
            switch(size(label,2))
                case 1
                    parameter = setfield(parameter,HEADLINE,{BF},label{1},...
                    str2num(tline(strfind(tline,'=')+1:end)));
                case 2
                    parameter = setfield(parameter,HEADLINE,{BF},label{1},label{2},...
                    str2num(tline(strfind(tline,'=')+1:end)));
                case 3
                    parameter = setfield(parameter,HEADLINE,{BF},label{1},label{2},label{3},...
                    str2num(tline(strfind(tline,'=')+1:end)));
                otherwise                  
            end
        end
    end
    
     % read till end of file
    if(~ischar(tline) & tline == -1)
%         disp(['Processed ' num2str(count) ' lines'])
        fclose(fid);
        break
    end

    count = count +1;
end
warning('on','MATLAB:dispatcher:InexactMatch')
% Please send suggestions for improvement to Martin Christian Hemmsen at 
% this email address: mah@elektro.dtu.dk
% Your contribution towards improving this code will be acknowledged in
% the "Changes" section.