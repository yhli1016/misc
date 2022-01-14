function lmn = dircos(dx)
%
% This function calculate the direction cosines used in Slater-Koster
% relation.
%
% Input:
%   dx: displacement vector pointing from atom1 to atom2, in Cartesian
%       coordinates
%
% Output:
%  lmn: direction cosines
%
    dx_mod = sqrt(sum(dx.^2));
    if (dx_mod <= 1.0e-9)
        lmn = 0.0 * dx;
    else
        lmn = dx / dx_mod; 
    end
