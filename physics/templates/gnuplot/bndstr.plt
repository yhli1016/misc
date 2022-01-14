## This script plots the common band structure.

reset

## set terminal
width_term = 6
height_term = 4
font_term = "font 'Arial, 12' fontscale 1.0"
outfile = "bndstr.pdf"

set term pdfcairo size width_term, height_term @font_term
set output outfile

## plot band structure
width_line = 2
width_border = 2
width_arrow = 2

# x-axis, location of high-symmetric k-points
x1 = 0.000000
x2 = 0.118334
x3 = 0.186654
x4 = 0.323295

# y-axis
ymin = -1.0
ymax = 2.5
dy = 0.5
nmytics = 5
efermi = -0.70560

# set box
set border lw width_border

# set axis
set xrange [x1:x4]
set xtics nomirror  ('{/Symbol G}' x1, 'M' x2, 'K' x3, '{/Symbol G}' x4)
unset x2tics

set yrange [ymin:ymax]
set ytics nomirror dy format "%4.1f"
set mytics nmytics
set ylabel 'Energy (eV)'
unset y2tics

# set legend
unset key

# draw auxiliary lines
set arrow from x2, ymin to x2, ymax nohead lw width_arrow
set arrow from x3, ymin to x3, ymax nohead lw width_arrow

# plot
plot for [i=2:161] 'bndstr.pd' u 1:(column(i)-efermi) w l lw width_line lc rgb 'blue'

## flush
set output
