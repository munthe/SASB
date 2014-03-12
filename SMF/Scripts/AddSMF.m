function [image] = SecondstageSMF( RFdata,resolution,SMFpath )
% Input
% RFdata, Firststage image
% resolution, [depth,number of lines]
% SMFpath, directory with SMF_line_[line number].mat files
% 
% Output
% image, Second stage image

image = zeros(resolution(1),92);
for line = 1:resolution(2)
    fprintf('Filtering %i line \n',line)
    load([SMFpath 'SMF_line_' num2str(line)], 'SMFline');
    imageLine = Secondstage_line(SMFline);
%     clear SMFline;
    image(:,line) = imageLine;
end
% image = imageLine;

function imageLine = Secondstage_line(filter)
    imageLine = zeros(resolution(1),1);
    for i = 1:92%resolution(1)
        point = RFdata( ...
                filter(i).index(1,1):filter(i).index(2,1) , ...
                filter(i).index(1,2):filter(i).index(2,2) );
        filtered = filter(i).filter .* point;
        imageLine(i) = sum(sum( filtered ));
    end
   
% SMF(coord).filter .* RFdata(SMF(coord).index(1,1):SMF(coord).index(2,1),SMF(coord).index(1,2):SMF(coord).index(2,2));
end

end

