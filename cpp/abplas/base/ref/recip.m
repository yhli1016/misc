function y = recip(x)
%
% x = [ a1x, a1y, a1z;
%       a2x, a2y, a2z;
%       a3x, a3y, a3z ];
%

    % NOTE: the volume here should have a sign, i.e. DO NOT
    % add any abs() here. Otherwise the results will be wrong
    % when vol < 0.
    a1  = x(1,:);
    a2  = x(2,:);
    a3  = x(3,:);
    vol = cross(a1, a2) * a3';
    b1  = 2 * pi/ vol * cross(a2, a3);
    b2  = 2 * pi/ vol * cross(a3, a1);
    b3  = 2 * pi/ vol * cross(a1, a2);
    y   = [b1; b2; b3];

    % Alternatively, you can use the algorithm below which
    % strictly follows the definition of reciprocal lattice
    % vectors.
    %y = zeros(3, 3);
    %product = 2 * pi * eye(3);
    %for ii = 1:3
    %    y(:,ii) = linsolve(x, product(:,ii));
    %end
    %y = y';
