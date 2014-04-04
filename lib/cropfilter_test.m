function [] = cropfilter_test()

SMFpath = '/data/cfudata3/mah/Spatial_matched_filter/SMF_3000x192_crop60dB/';
resolution = [3000,100];

for line = 1:resolution(2)
    tic
    fprintf('Cropping %i line ',line)
    load([SMFpath 'SMF_line_' num2str(line)], 'SMFline');
    SMFline = calcline(SMFline);
%     clear SMFline;
    save([SMFpath 'cropped/' 'SMF_line_' num2str(line)], 'SMFline');
    toc
end

function data = calcline(data)
    parfor i = 1:resolution(1)
        data(i).filter = cropfilter(data(i).filter,-20);
    end
end

end
