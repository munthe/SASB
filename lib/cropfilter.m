function SMF = cropfilter(SMF,dB)
% Sets filter values to zeros at the set threshold attenuation

SMFlog = logcompress(SMF);
SMF(SMFlog<dB) = 0;

end