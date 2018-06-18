## This script plots the transition matrix elements as the function of excited state index.

reset

## set terminal
width_term = 4
height_term = 3
font_term = "font 'Arial, 9' fontscale 1.0"
outfile = 'dip.pdf'

set term pdfcairo size width_term, height_term @font_term
set output outfile

## plot data
width_line = 2
width_border = 2
width_arrow = 2

# set box
set border lw width_border

# set axis
set xrange [-5:50]
set xtics 10 nomirror
set xlabel 'Excited state index'
unset x2tics

set yrange [0:400]
set ytics 100 nomirror
set mytics 2
set ylabel 'Dipole (a.u.)'
unset y2tics

# set legend
set key opaque reverse Left font 'Arial, 8'

# set histogram style
set style data histogram
set style histogram clustered gap 1
set style fill solid 0.4 border

# plot
plot "../0/15-abs_am/eigenvalues.dat" u 2 lc rgb 'red' t "E // AM", \
     "../0/15-abs_zz/eigenvalues.dat" u 2 lc rgb 'blue' t "E // ZZ", \
     "../0/15-abs_z/eigenvalues.dat"  u ($2*50) lc rgb 'dark-green' t "E // c"

## flush
set output
