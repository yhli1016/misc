function shiftlayer(c, d, z0, dr, layer1fn, layer2fn, outfn)
%
%        c: length of c-axis
%        d: desired interlayer distance
%       z0: desired geometric center along c-axis in fractional coordinates
%       dr: shift of layer_2 according to layer_1 in fractional coordinates
% layer1fn: file name of the data file containing the coordinates of layer1
% layer2fn: file name of the data file containing the coordinates of layer2
%    outfn: file name of the data file for output
%
% NOTE: layer1 is always assumed to be the bottom layer while layer2 is the
%       top layer.
%


    % initialization
    layer1 = load(layer1fn);
    layer2 = load(layer2fn);
    [ m1, n1 ] = size(layer1);
    [ m2, n2 ] = size(layer2);
    layer1new = layer1;
    layer2new = layer2;


    % shift layer2 such that the interlayer distance is d
    z1max = max(layer1(:,3));
    z2min = min(layer2(:,3));
    dz2 = z1max + d/c - z2min;
    layer2new = layer2new + [ 0, 0, dz2 ];


    % shift layer1 and layer2 such that the geometric center along c-axis is
    % at z0
    dz3 = z0 - ( sum(layer1new(:,3)) + sum(layer2new(:,3)) ) / ( m1 + m2 );
    layer1new = layer1new + [ 0, 0, dz3 ];
    layer2new = layer2new + [ 0, 0, dz3 ];


    % shift layer2 according to dr and check if any coordinate of layer2new
    % falls out of [0,1]
    layer2new = layer2new + dr;
    for i = 1:m2
        for j = 1:3
            if (layer2new(i,j) > 1.0)
                layer2new(i,j) = layer2new(i,j) - 1.0;
            elseif (layer2new(i,j) < 0.0)
                layer2new(i,j) = layer2new(i,j) + 1.0;
            end
        end
    end


    % write to file
    % FIRST LAYER 1
    outfid = fopen(outfn, 'w');
    for i = 1:m1
        fprintf(outfid, '%14.9f%14.9f%14.9f\n', layer1new(i,:));
    end
    % THEN LAYER 2
    for i = 1:m2
        fprintf(outfid, '%14.9f%14.9f%14.9f\n', layer2new(i,:));
    end
    fclose(outfid);
