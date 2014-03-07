function [ image ] = Second_Stage_SMF( RFdata, SMF )
%Second_Stage_SMF filters the first stage image
%   and saves it to the folder Second_Stage_RF_Data_Beamformed



% Saving does not work.


image = cellfun( @(X) ...
    sum(sum(X.*RFdata)),...
    SMF);

if(~isdir([savepath 'Second_Stage_RF_Data_Beamformed' filesep]))
            mkdir([savepath 'Second_Stage_RF_Data_Beamformed' filesep])
end
        RFdata = image;
        str = sprintf('%sSecond_Stage_RF_Data_Beamformed%s%sSecond_Stage_RF_Data_Beamformed.mat',savepath,filesep);
        save(str,'RFdata');  
end


% From BeamformationVer4 line 1073
% if(~isdir([savepath 'Second_Stage_RF_Data_Beamformed' filesep lower(type)]))
%            mkdir([savepath 'Second_Stage_RF_Data_Beamformed' filesep lower(type)])
%        end
%        RFdata = bf_image_bft3;
%        str = sprintf('%sSecond_Stage_RF_Data_Beamformed%s%s%sSecond_Stage_RF_Data_Beamformed.mat',savepath,filesep,lower(type),filesep);
%        save(str,'RFdata','useCaseParams','bft3','xmt');  