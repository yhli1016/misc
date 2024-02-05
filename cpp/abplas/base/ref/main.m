clear;

lat_vec = [1.0, 1.0, -1.0;
           2.0, -1.5, 1.0;
           0.0, 0.1, 2.1];
frac_coord = [0.5, 0.0, 0.0, 0.7, 6.0;
              1.0, 1.0, 0.0, 0.1, -0.1;
              0.0, 1.0, 1.0, -1.2, 0.7];
disp("lat_vec");
disp(lat_vec);
disp("frac_coord");
disp(frac_coord);

cart_coord = frac2cart(lat_vec', frac_coord')';
disp("cart_coord")
disp(cart_coord);

delta = cart2frac(lat_vec', cart_coord')' - frac_coord;
disp("abs(delta)");
disp(abs(delta));

recip_vec = recip(lat_vec')';
disp("recip_vec");
disp(recip_vec);
