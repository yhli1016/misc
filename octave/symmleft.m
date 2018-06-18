function y = symmleft(x)
%
% This function builds a whole symmetric function
% from the left part of x.
%
% x: N*2 array
%
    [m, n] = size(x);
    xleft = x(1:m-1,1);
    yleft = x(1:m-1,2);
    dx = x(2,1) - x(1,1);
    xright = xleft + (x(m,1) - x(1,1) + dx);
    yright = flipud(yleft);
    xfull = [x(:,1); xright];
    yfull = [x(:,2); yright];
    y = [xfull, yfull];
