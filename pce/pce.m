function y = pce(Eg, Ec)
    global AM15G;
    
    % determine the numerator
    [m ,n] = size(AM15G);
    imin = 1;
    for i = 1:m
        if (AM15G(i,1) >= Eg)
            imin = i;
            break;
        end
    end
    Numerator = trapz(AM15G(imin:m,1), AM15G(imin:m,2) ./  AM15G(imin:m,1));
    
    % determine the denorminator
    Denorminator = trapz(AM15G(:,1), AM15G(:,2));
    
    % calculate PCE
    y = 0.65 * (Eg - Ec - 0.3) * Numerator / Denorminator ;
