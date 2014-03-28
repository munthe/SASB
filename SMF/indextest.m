

load '/data/cfudata3/mah/Spatial_matched_filter/SMF_3000x192_crop10dB_line91-102/SMF_line_96/';

 l =(useCaseParams.scanparams(1).windowtissueq.y_tismax-useCaseParams.scanparams(1).windowtissueq.y_tismin);
 
 dep = 0.075;
 dep = dep -useCaseParams.scanparams(1).windowtissueq.y_tismin;
 tal = round((dep/l)*1000)
 
 filter1 = SMFline(tal).filter;
 figure(3)
 imagesc(filter1)
 
 RFdata1 = RFdata(SMFline(tal).index(1,1):SMFline(tal).index(2,1),SMFline(tal).index(1,2):SMFline(tal).index(2,2));
 figure(4)
 imagesc(RFdata1)
 
 test1 = RFdata1.*filter1;
 figure(5)
 imagesc(test1)