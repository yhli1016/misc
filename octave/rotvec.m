function y = rotvec(u, theta, x)
% Function rotating vectors in x about given axis u with angle theta.

% See https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
% for details.
% Note that the right-hand rule is used, i.e., a positive value means the angle
% is counter-clockwise.

% u: [ux, uy, uz], row vector defining the fixed axis
% theta: rotation angle in DEGREEs
% x = [v1x, v1y, v1z;
%      v2x, v2y, v2z;
%      ..............
%      vNx, vNy, vNz]

    % determine the rotation matrix
    ux = u(1); uy = u(2); uz = u(3);
    uprod = [0, -uz, uy; uz, 0, -ux; -uy, ux, 0];
    
    % The formula tensor_product(u, u) = u * u' in the URL is for column vector
    % while we have a row vector of u. So we have to transpose u here.
    utens = u' * u;
    theta2 = theta / 180 * pi;
    convmat = cos(theta2) * eye(3) + sin(theta2) * uprod + (1 - cos(theta2)) * utens;
    
    % rotate vectors
    [m, n] = size(x);
    y = zeros(m, n);
    for k = 1:m
        y(k,:) = (convmat * x(k,:)')';
    end
