function savemat(fnm, x, fmt)
    fid = fopen(fnm, 'w');
    [ m, n ] = size(x);
    for k = 1:m
        for l = 1:n
            fprintf(fid, fmt, x(k,l));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);