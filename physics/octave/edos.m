function y = edos(x, eng, sigma, flavor)
%
% This function calculates the density of excited states
% according to their energies.
%
%      y: calculated density of excited states, n*1 array
%      x: array containing energies to evaluate edos, n*1 array 
%    eng: energies of excited states, n*1 array
%  sigma: broadening parameter, float
% flavor: type of broadening funcition, 'gaussian' or 'lorentzian'
%         available
%
    y = zeros(size (x));
    N = length(eng);
    if (strcmp(flavor, 'gaussian'))
        for i = 1:N
            y = y + gaussian(x, eng(i), sigma);
        end
    elseif (strcmp(flavor, 'lorentzian'))
            for i = 1:N
            y = y + lorentzian(x, eng(i), sigma);
        end
    else
        disp('Unknown broadening type!');
        disp(flavor);
    end