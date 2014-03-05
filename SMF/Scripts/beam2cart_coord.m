function cart = mfr_beam2cart_coord (r, zy, zx)
%
% cartesian = beam2cart_coord (r, zy, zx)
%
% size(cartesian,1) = length(r)
% size(cartesian,2) = 3 [x y z]
%
% By MFR, Init version some time in August 2011
% Version 1.1 2011-11-05 MFR: Changed name from beam_coord_to_cart to beam2cart_coord
% Version 1.2 2011-11-07 MFR: Input var verification and changed NARGIN from 1 to 3. Removed the loop.
%

dimR = size(r);
dimZY = size(zy);
dimZX = size(zx);

if min(dimR) > 1,  error('r must be scalor or vector'); end
if min(dimZY) > 1, error('zy must be scalor or vector'); end
if min(dimZX) > 1, error('zx must be scalor or vector'); end

% force column vector
if dimR(1)  < dimR(2),   r  = r'; end;
if dimZY(1) < dimZY(2), zy = zy'; end;
if dimZX(1) < dimZX(2), zx = zx'; end;

if (max(dimR) ~= max(dimZY)) || (max(dimZY) ~= max(dimZX)), error('All inputs must have the same length.'); end

% do the calc
z = sqrt(r.^2./(1+ tan(zx).^2 + tan(zy).^2));
x = z.*tan(zx);
y = z.*tan(zy);

% output
cart = [x y z];

