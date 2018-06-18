## This script plots the absorption spectrum.

reset

## set terminal
width_term = 4
height_term = 3
font_term = "font 'Arial, 9' fontscale 1.0"
outfile = 'abs.pdf'

set term pdfcairo size width_term, height_term @font_term enhanced
set output outfile

## plot data
width_line = 2
width_border = 2
width_arrow = 2

# set key
set key reverse Left font 'Arial, 8'

# set box
set border lw width_border

# set axis
set xrange [2:6]
set xtics 1 nomirror
set mxtics 2
set xlabel 'Energy (eV)'
unset x2tics

set yrange [0:20]
set ytics 5 nomirror
set mytics 2
set ylabel '{/Symbol e}_2 (arb.unit)'
unset y2tics

# plot
plot "../0/15-abs_am/absorption_eh.dat" u 1:2 w l lw width_line lc rgb 'red' t "E // AM", \
     "../0/15-abs_zz/absorption_eh.dat" u 1:2 w l lw width_line lc rgb 'blue' t "E // ZZ", \
     "../0/15-abs_z/absorption_eh.dat"  u 1:($2*50) w l lw width_line lc rgb 'dark-green' t "E // c"

## flush
set output
