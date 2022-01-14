clear;

% total simulation time in SECONDS
duration = 1.0 * 10^-12;
datfile = "max.dat";

dat = load(datfile);
eng = dat(:,3);
[eng_fft, freq] = fft_wrapper(eng, duration);

c = 299792458;
wavenumber = freq / c / 100;
plot(wavenumber, abs(eng_fft));
set(gca, 'xlim', [100,2000]);
set(gca, 'ylim', [0,60]);
