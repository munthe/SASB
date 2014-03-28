function [] = cropfilter_test()

SMFpath = '/usr/local/s103303/SMF_1000x100/';
resolution = [1000,100];

for line = 1:resolution(2)
    fprintf('Cropping %i line \n',line)
    load([SMFpath 'SMF_line_' num2str(line)], 'SMFline');
    SMFline = calcline(SMFline);
%     clear SMFline;
    save([SMFpath 'cropped/' 'SMF_line_' num2str(line)], 'SMFline');
end

function data = calcline(data)
    for i = 1:resolution(1)
        data(i).filter = cropfilter(data(i).filter,-60);
    end
end

end
