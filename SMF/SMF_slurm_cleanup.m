% 
% 

%% Load filter lines
SMF = cell(resolution(1),length(par.line));
for line = par.line
    load([savepath tmpdir 'SMF_line_' num2str(line)]);
    SMF(:,line) = SMFline;
end

%% Save filter
save([savepath 'SMF'],'SMF')

%% Clean up
rmdir([savepath tmpdir],'s')
mkdir([savepath tmpdir])

