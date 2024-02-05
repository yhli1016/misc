function y = cart2frac(latvec, x)
%
% latvec = [ a1x, a1y, a1z;
%            a2x, a2y, a2z;
%            a3x, a3y, a3z ]
%
% x = [ atom1x, atom1y, atom1z;
%       atom2x, atom2y, atom2z;
%       .......................
%       atomNx, atomNy, atomNz ]
%
    [ m, n ] = size(x);
    y = zeros(m, n);
    convmat = inv(latvec');
    for k = 1:m
        y(k,:) = (convmat * x(k,:)')';
    end
