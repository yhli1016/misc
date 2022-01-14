function y = eos(vol, eng, N)
%
% This function fits volume-energy relation and find
% the equilibrium volume.
%
    p = polyfit(vol, eng, N);
    dp = polyder(p);
    y = roots(dp);
    xfi = linspace(min(vol), max(vol), 100);
    yfi = polyval(p, xfi);
    plot(vol, eng, "+");
    hold on;
    plot(xfi, yfi, "-1");
