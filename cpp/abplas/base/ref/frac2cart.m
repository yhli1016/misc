function y = frac2cart(latvec, x)
%
% latvec = [ a1x, a1y, a1z;
%            a2x, a2y, a2z;
%            a3x, a3y, a3z ]
%
% x = [ atom11, atom12, atom13;
%       atom21, atom22, atom23;
%       .......................
%       atomN1, atomN2, atomN3 ]
%
    [ m, n ] = size(x);
    y = zeros(m, n);
    convmat = latvec';
    for k = 1:m
        y(k,:) = (convmat * x(k,:)')';
    end
