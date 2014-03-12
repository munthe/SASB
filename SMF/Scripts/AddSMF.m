function [RFdata2] = AddSMF( RFdata,resolution,SMFpath )
% Adds a Spatial Matched Filter to a First Stage image
%
% Input
% RFdata, Firststage image
% resolution, [depth,number of lines]
% SMFpath, directory with SMF_line_[line number].mat files
% 
% Output
% image, Second stage image

RFdata2 = zeros(resolution);
for line = 1:resolution(2)
    fprintf('Filtering %i line \n',line)
    load([SMFpath 'SMF_line_' num2str(line)], 'SMFline');
    scanline = Secondstage_line(SMFline);
%     clear SMFline;
    RFdata2(:,line) = scanline;
end
% image = scanline;

function scanline = Secondstage_line(filter)
    scanline = zeros(resolution(1),1);
    for i = 1:resolution(1)
        f = crop2view(RFdata,filter(i));
        point = RFdata( ...
                f.index(1,1):f.index(2,1) , ...
                f.index(1,2):f.index(2,2) );
        filtered = f.filter .* point;
        scanline(i) = sum(sum( filtered ));
    end
   
end

end

function filter = crop2view(RFdata,filter)
% Checks if the filter index goes beyond the depth of the RFdata, if so
% crop the filter and update the filter index.
delta = filter.index(2,1) - size(RFdata,1);
if delta > 0
    filter.index(2,1) = size(RFdata,1);
    filter.filter = filter.filter(1:end-delta,:);
end
end
