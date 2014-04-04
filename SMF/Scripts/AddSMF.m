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

RFdata2 = zeros(resolution(2,:));
for line = resolution(1,2):resolution(2,2)/2
    fprintf('Filtering line %i and %i ... ',line,resolution(2,2)-line+1)
    tic
    load([SMFpath 'SMF_line_' num2str(line)], 'SMFline');
    scanline_l = Secondstage_line_l(SMFline);
%     clear SMFline;
    RFdata2(:,line) = scanline_l;
%     SMFline(line).filter = arrayfun(@(l)(fliplr(SMFline(l).filter)),1:resolution(2,1),'UniformOutput',false);
%     SMFline(line).index(:,2) = [resolution(2,2)-SMFline(line).index(2,2)+1;resolution(2,2)-SMFline(line).index(2,1)+1]; 
    scanline_r = Secondstage_line_r(SMFline);
    RFdata2(:,end-line+1) = scanline_r;
    fprintf('took %f seconds \n',toc)
end
% image = scanline;

function scanline = Secondstage_line_l(filter)
    scanline = zeros(resolution(2,1),1);
    parfor i = resolution(1,1):resolution(2,1)
        f = crop2view(RFdata,filter(i));
        point = RFdata( ...
                f.index(1,1):f.index(2,1) , ...
                f.index(1,2):f.index(2,2) );
        filtered = f.filter .* point;
        scanline(i) = sum(sum( filtered ));
    end   
end
function scanline = Secondstage_line_r(filter)
    scanline = zeros(resolution(2,1),1);
    parfor i = resolution(1,1):resolution(2,1)
        filter(i).filter = fliplr(filter(i).filter);
        filter(i).index(:,2) = [resolution(2,2)-filter(i).index(2,2)+1;resolution(2,2)-filter(i).index(1,2)+1];
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

