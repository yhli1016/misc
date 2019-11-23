clear;

% Lower and upper bound of band gap of donor and conduction bands offset
Eg_min = 0.5;
Eg_max = 3.5;
Ec_min = 0.0;
Ec_max = 1.0;

% load solar spectrum and calculate the denorminator
global AM15G;
AM15G = load('AM15G.DAT');

% evaluate ita on a coarse grid
NGrid = 25;
Eg = linspace(Eg_min, Eg_max, NGrid);
Ec = linspace(Ec_min, Ec_max, NGrid);
ita = zeros(NGrid * NGrid, 3);
for i = 1:NGrid
    for j = 1:NGrid
        k = (i - 1) * NGrid + j;
        ita(k,1) = Eg(i);
        ita(k,2) = Ec(j);
        ita(k,3) = pce(Eg(i), Ec(j));
    end
end

% interpolate onto a finer grid and plot
Eg_fi = linspace(Eg_min, Eg_max, 2 * NGrid);
Ec_fi = linspace(Ec_min, Ec_max, 2 * NGrid);
[XX, YY] = meshgrid(Eg_fi, Ec_fi);
ZZ = griddata(ita(:,1), ita(:,2), ita(:,3), XX, YY);
contourf(XX, YY, ZZ);
colorbar;
