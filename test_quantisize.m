% Script to test quantize funtion

han=hanning(192);
levels = 16;

% Levels in powers of two
expo = [0,2.^[0:levels-2]];
expo = expo/max(max(expo));

han_expo = quantize(han,expo);

% Linear levels 
han_lin = quantize(han,0:1/levels:1);

figure(1)
hold on
plot(han_expo,'r')
plot(han_lin,'b')