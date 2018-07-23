clear;

prefix = "eps1_8.1_10";

# Load data
dat = load([prefix, "/eps2.dat"]);
xco = dat(:, 1);
yco = dat(:, 3);
zco = dat(:, 5);

# Interpolate to a denser grid
xfi = linspace(min(xco), max(xco), 100);
yfi = linspace(min(yco), max(yco), 100);
[xx, yy] = meshgrid(xfi, yfi);
zz = griddata(xco, yco, zco, xx, yy);

# Plot
graphics_toolkit("gnuplot");
contourf(xx, yy, zz);
colormap("jet");
colorbar();
set(gca, "xlabel", "Thickness (nm)");
set(gca, "ylabel", "Dielectric constant");
set(gca, "title", "Binding energy");
print(gcf, [prefix, "/eb2_oct.pdf"]);