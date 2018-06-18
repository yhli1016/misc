function y = lorentzian(x, mu, sigma)
    part_a = 1.0 / (pi * sigma);
    part_b = sigma ^ 2 ./ ((x - mu) .^ 2 + sigma ^ 2);
    y = part_a * part_b;