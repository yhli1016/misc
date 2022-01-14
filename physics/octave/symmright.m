function y = symmright(x)
%
% This function builds a whole symmetric function
% from the right part of x.
%
% x: N*2 array
%
    [m, n] = size(x);
    xright = x(2:m,1);
    yright = x(2:m,2);
    dx = x(2,1) - x(1,1);
    xleft = xright - (x(m,1) - x(1,1) + dx);
    yleft = flipud(yright);
    xfull = [xleft; x(:,1)];
    yfull = [yleft; x(:,2)];
    y = [xfull, yfull];
