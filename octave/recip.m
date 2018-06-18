function y = recip(x)
%
% x = [ a1x, a1y, a1z;
%       a2x, a2y, a2z;
%       a3x, a3y, a3z ];
%
    a1  = x(1,:);
    a2  = x(2,:);
    a3  = x(3,:);
    vol = abs(cross(a1, a2) * a3');
    b1  = 2*pi/vol * cross(a2, a3);
    b2  = 2*pi/vol * cross(a3, a1);
    b3  = 2*pi/vol * cross(a1, a2);
    y   = [ b1; b2; b3 ];
