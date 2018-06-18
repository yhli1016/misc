clear;

duration = 10;
NSample = 5000;
t = linspace(0, duration, NSample)';
yt = 0.1 + 0.5 * sin(2*pi*1*t) + sin(2*pi*2*t) + 0.5 * sin(2*pi*3*t) + 0.1 * sin(2*pi*4*t);
[yf, f] = fft_wrapper(yt, duration);

plot(f, abs(yf));
set(gca, 'xlim', [0,10]);
set(gca, 'xtick', [0:1:10]);