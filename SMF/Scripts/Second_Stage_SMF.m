function [ image ] = Second_Stage_SMF( RFdata, SMF, useCaseParams )
%Second_Stage_SMF filters the first stage image
%   and saves it to the folder Second_Stage_RF_Data_Beamformed

resolution = [100 , 100];
SMFpath = 'tmp_SMF_100x100/';
image = AddSMF(RFdata,resolution,SMFpath);

savepath = '/home/s113291/SASB/SMF/';
%savepath = loadpath;

% remove any old directory
        if(isdir([savepath 'Second_Stage_RF_Data_Beamformed' ]))
           rmdir([savepath 'Second_Stage_RF_Data_Beamformed'],'s')
        end
        
if(~isdir([savepath 'Second_Stage_RF_Data_Beamformed']))
           mkdir([savepath 'Second_Stage_RF_Data_Beamformed'])
end

      RFdata = image;
      str = sprintf('%sSecond_Stage_RF_Data_Beamformed',savepath, filesep);
      save(str,'RFdata','useCaseParams');  
