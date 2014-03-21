function image = logcompress(RFdata)

image = abs(hilbert(RFdata));
image = 20*log10(image/max(max(image)));

end
