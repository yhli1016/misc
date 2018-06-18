## This script plots the absorption spectra for many-systems in one figure.

reset

## set terminal
width_term = 10
height_term = 6
font_term = "font 'Arial, 8' fontscale 1.0"
outfile = 'abs.pdf'

set term pdfcairo size width_term, height_term @font_term
set output outfile
set multiplot layout 3,3 columnsfirst

## plot data
width_line = 2
width_border = 2
width_arrow = 2

# set box
set border lw width_border

# set x-axis
set xrange [2:6]
set xtics 1
set mxtics 2
set xlabel 'Energy (eV)'

# set y-axis and plot
set mytics 2
set ylabel '{/Symbol e}_2 (arb. unit)'

ymax = 20
dy = 5
set yrange [0:ymax]
set ytics dy
set arrow from 2.83664,0 to 2.83664,ymax nohead lw width_arrow lc rgb 'red' dt 4
set arrow from 2.76599,0 to 2.76599,ymax nohead lw width_arrow lc rgb 'blue' dt 4
do for [pref in "am-10.0 am-5.5 am+5.0"] {
    plot pref."/15-abs_am/absorption_eh.dat" u 1:2 w l lw width_line lc rgb 'red' t "AM", \
         pref."/15-abs_zz/absorption_eh.dat" u 1:2 w l lw width_line lc rgb 'blue' t "ZZ"
}

ymax = 10
dy = 2.5
set yrange [0:ymax]
set ytics dy
unset arrow
set arrow from 2.83664,0 to 2.83664,ymax nohead lw width_arrow lc rgb 'red' dt 4
set arrow from 2.76599,0 to 2.76599,ymax nohead lw width_arrow lc rgb 'blue' dt 4
do for [pref in "zz-10.0 zz-5.5"] {
    plot pref."/15-abs_am/absorption_eh.dat" u 1:2 w l lw width_line lc rgb 'red' t "AM", \
         pref."/15-abs_zz/absorption_eh.dat" u 1:2 w l lw width_line lc rgb 'blue' t "ZZ"
}

ymax = 20
dy = 5
set yrange [0:ymax]
set ytics dy
unset arrow
set arrow from 2.83664,0 to 2.83664,ymax nohead lw width_arrow lc rgb 'red' dt 4
set arrow from 2.76599,0 to 2.76599,ymax nohead lw width_arrow lc rgb 'blue' dt 4
do for [pref in "zz+5.0"] {
    plot pref."/15-abs_am/absorption_eh.dat" u 1:2 w l lw width_line lc rgb 'red' t "AM", \
         pref."/15-abs_zz/absorption_eh.dat" u 1:2 w l lw width_line lc rgb 'blue' t "ZZ"
}

ymax = 15
dy = 5
set yrange [0:ymax]
set ytics dy
unset arrow
set arrow from 2.83664,0 to 2.83664,ymax nohead lw width_arrow lc rgb 'red' dt 4
set arrow from 2.76599,0 to 2.76599,ymax nohead lw width_arrow lc rgb 'blue' dt 4
do for [pref in "bi-6.5 bi-3.5 bi+5.0"] {
    plot pref."/15-abs_am/absorption_eh.dat" u 1:2 w l lw width_line lc rgb 'red' t "AM", \
         pref."/15-abs_zz/absorption_eh.dat" u 1:2 w l lw width_line lc rgb 'blue' t "ZZ"
}

## flush
unset multiplot
set output
