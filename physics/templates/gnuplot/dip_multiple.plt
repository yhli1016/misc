## This script plots the dipole for many-systems in one figure.

reset

## set terminal
width_term = 10
height_term = 6
font_term = "font 'Arial, 8' fontscale 1.0"
outfile = 'dip.pdf'

set term pdfcairo size width_term, height_term @font_term
set output outfile
set multiplot layout 3,3 columnsfirst

## plot data
width_line = 2
width_border = 2
width_arrow = 2

# set box
set border lw width_border

# set legend
set key opaque

# set x-axis
set xrange [-5:50]
set xtics 10 nomirror
set xlabel 'Excited state index'
unset x2tics

# set histogram style
set style data histogram
set style histogram clustered gap 1
set style fill solid 0.4 border

# set y-axis and plot
set mytics 2
set ylabel 'Dipole (a.u.)'

ymax = 400
dy = 100
set yrange [0:ymax]
set ytics dy nomirror
unset y2tics
do for [pref in "am-10.0 am-5.5"] {
    plot pref."/15-abs_am/eigenvalues.dat" u 2  lc rgb 'red' t "AM", \
         pref."/15-abs_zz/eigenvalues.dat" u 2  lc rgb 'blue' t "ZZ"

}

ymax = 500
dy = 100
set yrange [0:ymax]
set ytics dy
do for [pref in "am+5.0"] {
    plot pref."/15-abs_am/eigenvalues.dat" u 2  lc rgb 'red' t "AM", \
         pref."/15-abs_zz/eigenvalues.dat" u 2  lc rgb 'blue' t "ZZ"

}


ymax = 100
dy = 20
set yrange [0:ymax]
set ytics dy
do for [pref in "zz-10.0 zz-5.5"] {
    plot pref."/15-abs_am/eigenvalues.dat" u 2  lc rgb 'red' t "AM", \
         pref."/15-abs_zz/eigenvalues.dat" u 2  lc rgb 'blue' t "ZZ"

}


ymax = 500
dy = 100
set yrange [0:ymax]
set ytics dy
do for [pref in "zz+5.0"] {
    plot pref."/15-abs_am/eigenvalues.dat" u 2  lc rgb 'red' t "AM", \
         pref."/15-abs_zz/eigenvalues.dat" u 2  lc rgb 'blue' t "ZZ"

}

ymax = 400
dy = 100
set yrange [0:ymax]
set ytics dy
do for [pref in "bi-6.5 bi-3.5"] {
    plot pref."/15-abs_am/eigenvalues.dat" u 2  lc rgb 'red' t "AM", \
         pref."/15-abs_zz/eigenvalues.dat" u 2  lc rgb 'blue' t "ZZ"

}

ymax = 500
dy = 100
set yrange [0:ymax]
set ytics dy
do for [pref in "am+5.0"] {
    plot pref."/15-abs_am/eigenvalues.dat" u 2  lc rgb 'red' t "AM", \
         pref."/15-abs_zz/eigenvalues.dat" u 2  lc rgb 'blue' t "ZZ"

}

## flush
unset multiplot
set output
