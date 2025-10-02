set terminal pngcairo size 1200,400 enhanced font 'Verdana,21'
set output 'fe_result.png'
set xlabel 'x'
set ylabel 'u(x)'
set grid
plot 'fe_result_fem.dat' using 1:2 with linespoints lt 1 lw 2 pt 7 ps 1.2 title 'FEM', \
     'fe_result_analytic.dat' using 1:2 with linespoints lt 2 lw 2 pt 3 ps 1.2 title 'Analytical'