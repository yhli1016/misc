clear;

lat_vec = [1.0, 1.0, -1.0;
           2.0, -1.5, 1.0;
           0.0, 0.1, 2.1];
frac_coord = [0.5, 0.0, 0.0, 0.7, 6.0;
              1.0, 1.0, 0.0, 0.1, -0.1;
              0.0, 1.0, 1.0, -1.2, 0.7];
disp("\nlat_vec");
disp(transpose(lat_vec));
disp("\nfrac_coord");
disp(transpose(frac_coord));

cart_coord = zeros(3, 5);
for i = 1:5
    cart_coord(:, i) = lat_vec * frac_coord(:,i);
end
disp("\ncart_coord")
disp(transpose(cart_coord));

conv_mat = inv(lat_vec);
for i = 1:5
    frac_coord(:, i) = conv_mat * cart_coord(:,i);
end
disp("\nfrac_coord");
disp(transpose(frac_coord));