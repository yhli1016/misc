function y = normcoord(x)
% Function checks for fractional coordinates out of [0, 1]
    y = x;
    [m, n] = size(x);
    for i = 1:m
        for j = 1:n
            while (y(i,j) < 0.0)
                y(i,j) = y(i,j) + 1.0;
            end
            while (y(i,j) > 1.0)
                y(i,j) = y(i,j) - 1.0;
            end
        end
    end
