## This script plots the absorption spectra for many-systems in one figure.

reset

## set terminal
width_term = 8
height_term = 6
font_term = "font 'Arial, 8' fontscale 1.0"
outfile = 'bi.pdf'

set term pdfcairo size width_term, height_term @font_term
set output outfile
set multiplot layout 3,2 columnsfirst

## plot data
width_line = 2
width_border = 2
width_arrow = 2

# set box
set border lw width_border

# set key
set key Left reverse opaque

## First the absorption spectra
#
# set x-axis
set xrange [1:6]
set xtics 1 nomirror
set mxtics 2
set xlabel 'Energy (eV)'
unset x2tics

# set y-axis
ymax = 15
dy = 5

set yrange [0:ymax]
set ytics dy nomirror
set mytics 2
set ylabel '{/Symbol e}_2 (arb. unit)'
unset y2tics

set arrow from 2.84,0 to 2.84,ymax nohead lw width_arrow lc rgb 'red' dt 4
set arrow from 2.77,0 to 2.77,ymax nohead lw width_arrow lc rgb 'blue' dt 4
set arrow from 2.32,0 to 2.32,ymax nohead lw width_arrow lc rgb 'dark-green' dt 4

# plot
do for [pref in "bi-6.5 bi-3.5 bi+5.0"] {
    plot "../".pref."/15-abs_am/absorption_eh.dat" u 1:2 w l lw width_line lc rgb 'red' t "E // AM",\
         "../".pref."/15-abs_zz/absorption_eh.dat" u 1:2 w l lw width_line lc rgb 'blue' t "E // ZZ",\
         "../".pref."/15-abs_z/absorption_eh.dat"  u 1:($2*50) w l lw width_line lc rgb 'dark-green' t "E // c"
}


## Then the dipole moments
#
# overwrite previous settings
unset arrow

# set histogram style
set style data histogram
set style histogram clustered gap 1
set style fill solid 0.4 border

# set x-axis
set xrange [-5:50]
set xtics 10 nomirror
set xlabel 'Excited state index'
unset x2tics

# set y-axis
ymax = 300
dy = 100
pref = "bi-6.5"

set yrange [0:ymax]
set ytics dy nomirror
set mytics 2
set ylabel 'Dipole (a.u.)'
unset y2tics

# plot
plot "../".pref."/15-abs_am/eigenvalues.dat" u 2 lc rgb 'red' t "E // AM",\
     "../".pref."/15-abs_zz/eigenvalues.dat" u 2 lc rgb 'blue' t "E // ZZ",\
     "../".pref."/15-abs_z/eigenvalues.dat" u ($2*50) lc rgb 'dark-green' t "E // c"

ymax = 400
pref = "bi-3.5"
set yrange [0:ymax]
plot "../".pref."/15-abs_am/eigenvalues.dat" u 2 lc rgb 'red' t "E // AM",\
     "../".pref."/15-abs_zz/eigenvalues.dat" u 2 lc rgb 'blue' t "E // ZZ",\
     "../".pref."/15-abs_z/eigenvalues.dat" u ($2*50) lc rgb 'dark-green' t "E // c"

ymax = 500
pref = "bi+5.0"
set yrange [0:ymax]
plot "../".pref."/15-abs_am/eigenvalues.dat" u 2 lc rgb 'red' t "E // AM",\
     "../".pref."/15-abs_zz/eigenvalues.dat" u 2 lc rgb 'blue' t "E // ZZ",\
     "../".pref."/15-abs_z/eigenvalues.dat" u ($2*50) lc rgb 'dark-green' t "E // c"

## flush
unset multiplot
set output
