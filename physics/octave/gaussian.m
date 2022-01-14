function y = gaussian(x, mu, sigma)
    part_a = 1.0 / (sigma * sqrt(2 * pi));
    part_b = exp(-(x - mu) .^ 2 / (2 * sigma ^ 2));
    y = part_a * part_b;