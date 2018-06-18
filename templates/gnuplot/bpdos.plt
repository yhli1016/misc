## This script plots band structure and PDOS in one figure.

reset

## set terminal
width_term  = 6.6
height_term = 4.4
font_term = "font 'Arial, 10' fontscale 1.0"
outfile = 'bpdos.pdf'

set term pdfcairo size width_term, height_term @font_term
set output outfile
set multiplot

## plot band structure
width_line   = 2
width_border = 2
width_arrow  = 2
margin_bottom = 3

x1 = 0.000000
x2 = 0.619194
x3 = 0.979243
x4 = 1.482994
x5 = 1.843043

ymin = -15
ymax = 5
dy = 5
nmytics = 5

# set size, origin and margins
set size 0.6, 1.0
set origin 0.0, 0.0
set bmargin margin_bottom

# set box
set border lw width_border

# set axis
set xrange [x1:x5]
set xtics nomirror  ('{/Symbol G}' x1, 'M' x2, 'X' x3, '{/Symbol G}' x4, 'Y' x5)
unset x2tics

set yrange [ymin:ymax]
set ytics nomirror dy out
set mytics nmytics
set ylabel 'Energy (eV)'
unset y2tics

# set legend
unset key

# draw auxiliary lines
set arrow from x2, ymin to x2, ymax nohead lw width_arrow
set arrow from x3, ymin to x3, ymax nohead lw width_arrow
set arrow from x4, ymin to x4, ymax nohead lw width_arrow

# plot
plot for [i=2:61] 'bndstr.pd' u 1:i w l lw width_line lc rgb 'blue'

## plot PDOS
pdos_min = 0
pdos_max = 6
dpdos = 2
nmxtics = 2

reset

# set size, origin and margins
set size 0.4, 1.0
set origin 0.6, 0.0
set bmargin margin_bottom

# set box
set border lw width_border

# set axis
set xrange [pdos_min:pdos_max]
set xtics out nomirror dpdos
set mxtics nmxtics
unset x2tics

set yrange [ymin:ymax]
set ytics nomirror dy out
set mytics nmytics
set xlabel "PDOS (1/eV)" offset char 0, 0.5
unset y2tics

# set legend
set key right top reverse Left samplen 2 opaque

# plot
plot 'P3s.pd' u 2:1 w l lw width_line lc rgb 'red' t 'P 3s',\
     'P3p.pd' u 2:1 w l lw width_line lc rgb 'dark-green' t 'P 3p',\
     'P3d.pd' u 2:1 w l lw width_line lc rgb 'blue' t 'P 3d'

## flush
unset multiplot
set output
