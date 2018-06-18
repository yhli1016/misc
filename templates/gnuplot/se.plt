## This script plots the projection-resolved band structure.

reset

## set terminal
width_term = 5
height_term = 5
font_term = "font 'Arial, 14' fontscale 1.0"
outfile = "se.pdf"

set term pdfcairo size width_term, height_term @font_term
set output outfile

## plot band structure
width_line = 2
width_border = 2
width_arrow = 2

# x-axis, location of high-symmetric k-points
x1 = 0.000000
x2 = 0.473346
x3 = 0.746630
x4 = 1.293203

# y-axis
ymin = -6
ymax = 6
dy = 2
nmytics = 4

# energy shift for fermi level
efermi = -0.81740

# set box
set border lw width_border

# set axis
set xrange [x1:x4]
set xtics nomirror  ('{/Symbol G}' x1, 'M' x2, 'K' x3, '{/Symbol G}' x4)
unset x2tics

set yrange [ymin:ymax]
set ytics nomirror dy
set mytics nmytics
set ylabel 'Energy (eV)'
unset y2tics

# set legend
unset key

# draw auxiliary lines
set arrow from x2, ymin to x2, ymax nohead lw width_arrow
set arrow from x3, ymin to x3, ymax nohead lw width_arrow

# plot
plot 'se.pz.dat' u 1:($2-efermi) w l lw width_line lc rgb 'gray', \
'se.pz.dat' u 1:($2-efermi):3 w p lc rgb 'red' ps var pt 7, \
'se.pxy.dat' u 1:($2-efermi):($3*0.6) w p lc rgb 'blue' ps var pt 7

## flush
set output
