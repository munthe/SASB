function [ image ] = Second_Stage_SMF( RFdata, SMFpath,resolution,linespan,useCaseParams )
%Second_Stage_SMF filters the first stage image
%   and saves it to the folder Second_Stage_RF_Data_Beamformed

image = AddSMF(RFdata,[1,linespan(1);resolution(1),linespan(2)],SMFpath);

savepath = '../Figures_SMF';
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
