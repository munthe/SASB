function tgc = cfu_tgc_curve (size_data, dr, f0, loss_cm)
%
% returns a TGC-curve
%
% tgc = cfu_tgc_curve (size_data, dr, f0 [,loss_cm = 0.5 dB/(MHz*cm)])
%
% Input:
%   size-data:  Vector with the size of measured data in each dim. Ex.: size(samples) 
%   dr:         distance between sample points in meter
%   f0:         center frequency in Hz
%   loss_cm:    One-way signal loss in dB/(cm*MHz). Standard: 0.5dB/(MHz*cm)
%
% Output:
%   tgc_curve.dB:  TGC-curve in logarithmic (dB) scale
%   tgc_curve.lin: TGC-curve in linear scale
%
% Example: 
%    tgc = tgc_curve (size(samples), sys.dr, sys.f0)
%
% By Morten Fischer Rasmussen
% Version 1.0, 2011-12-02, Init version.
% Version 1.1, 2011-12-02, Renamed to tgc_curve. 
% Version 1.2, 2013-09-19, Renamed to cfu_tgc_curve. 
%


if nargin < 4
    loss_cm = 0.5; %[dB/(MHz*cm)]
end

no_samples = size_data(1);
tgc.dB = zeros(size_data);

loss_m = 2*loss_cm*100*f0/1e6;
spatial = (0:no_samples-1)'*dr;

tgc.dB  = repmat(spatial*loss_m, [1 size_data(2:end)]);
tgc.lin = 10.^(tgc.dB/20);



