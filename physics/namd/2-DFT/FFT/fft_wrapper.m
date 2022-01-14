function [y, f] = fft_wrapper(x, duration)
% Input:
%   x: data sampled in time-domain
%   duration: duration of sampling
%
% Output:
%   y: transformed data in ORDINARY frequency-domain
%   f: ORDINARY frequency, NOT the angular frequency omega!

    ndata = length(x);
    f = linspace(0, ndata-1, ndata)' / duration;
    
    % direct Fourier transform
    %y = zeros(ndata, 1);
    %t = linspace(0, ndata-1, ndata)' / ndata * duration;
    %for k = 1:ndata
    %    v = exp(-i * 2*pi * f(k) * t); 
    %    y(k) = v' * x;
    %end

    % built-in FFT
    y = fft(x);
