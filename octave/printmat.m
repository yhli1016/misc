function printmat(x, fmt)
    [ m, n ] = size(x);
    for k = 1:m
        for l = 1:n
            fprintf(fmt, x(k,l));
        end
        fprintf('\n');
    end
