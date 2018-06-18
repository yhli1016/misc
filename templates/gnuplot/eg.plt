## This script plots the Eg-strain relation.

reset

## set terminal
width_term = 8
height_term = 8
font_term = "font 'Arial, 14' fontscale 1.0"
outfile = 'eg.pdf'

set term pdfcairo size width_term, height_term @font_term
set output outfile

# plot Eg-strain relation
width_line = 2
width_border = 2
width_arrow = 2

# set box
set border lw width_border

# set axis
set xrange [-10.5:10.5]
set xtics nomirror 2
set mxtics 4
set xlabel 'Strain (%)'
unset x2tics

set ytics nomirror
set mytics 5
set ylabel 'E_g (eV)'
unset y2tics

# set legend
set key Left reverse nobox

# plot
plot 'eg.dat' u 1:2 w lp pt 6 lc rgb 'red'        t 'AM-Indirect', 'eg.dat' u 1:3 w lp pt 4 lc rgb 'red'        t 'AM-Direct',\
     'eg.dat' u 1:4 w lp pt 6 lc rgb 'blue'       t 'ZZ-Indirect', 'eg.dat' u 1:5 w lp pt 4 lc rgb 'blue'       t 'ZZ-Direct',\
     'eg.dat' u 1:6 w lp pt 6 lc rgb 'dark-green' t 'BI-Indirect', 'eg.dat' u 1:7 w lp pt 4 lc rgb 'dark-green' t 'BI-Direct'

## flush
set output
