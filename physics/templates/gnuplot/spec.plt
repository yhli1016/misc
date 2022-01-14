## This script plots the energy diagram of atoms.

reset

## set terminal
width_inch = 4
height_inch = 3
font_file = "font 'Arial, 10'"
outfile = 'spec.pdf'

set term pdfcairo size width_inch, height_inch @font_file
set output outfile

## plot data
width_line = 2
width_border = 2
width_arrow = 2

# set box
set border 2 lw  width_border

# set axis
set xrange [-1:6]
unset xtics

set yrange [-25:-5]
set ytics nomirror out
set mytics 5
set ylabel 'Energy (eV)'
unset y2tics

# set legend
unset key

# define line styles
set style line 1 lw width_line lc 'red'
set style line 2 lw width_line lc 'blue'

# plot
plot 'dia.dat' u 1:2 w l ls 1, 'dia.dat' u 3:4 w l ls 2

## flush output
set output
