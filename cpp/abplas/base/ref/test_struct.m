clear;

lat_vec = [
  1.2333 3.4444 1.0;
  1.5   2.7 2.1;
  1.0   3.9 1.2];

frac_coord = [
    1.0 0.5 0.0;
    2.0 1.0 0.3;
    0.1 1.2 0.1];

cart_coord = frac2cart(lat_vec, frac_coord);
disp(cart_coord);
