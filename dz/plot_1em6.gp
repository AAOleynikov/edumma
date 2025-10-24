# plot_phi_selected.gp
# data.dat:  idx t phi1 phi3 phi4
# res.tb9:   Time phi1 phi3 phi4

datafile1 = "1em6.dat"
datafile2 = "1em6.tb9"

img_w = 1200
img_h = 800

set datafile separator whitespace
set format x "%g"
set format y "%g"
set xlabel "t (s)"
set ylabel "phi (units)"
set grid
set style line 1 lw 2 pointsize 0.5
set style line 2 lw 2 dashtype 2 pointsize 0.5

# Массив колонок: phi1 -> col 3, phi3 -> col 4, phi4 -> col 5 для data.dat
# А для res.tb9 — Time = 1, phi(i) начиная с 2
names = "phi1 phi3 phi4"
cols_dat = "3 4 5"
cols_tb9 = "2 3 4"

do for [i=1:3] {
    col_dat = int(word(cols_dat, i))
    col_tb9 = int(word(cols_tb9, i))
    nm = word(names, i)
    outfile = sprintf("%s_1em6.png", nm)
    titletext = sprintf("%s", nm)

    set terminal pngcairo size img_w,img_h enhanced font "Arial,12"
    set output outfile

    set title titletext

    plot \
        datafile2 using 1:col_tb9 with linespoints ls 2 title sprintf("%s (ПА9)", nm), \
        datafile1 using 2:col_dat with linespoints ls 1 title sprintf("%s (Написанная программа)", nm)
        

    set output
}

# После выполнения появятся phi1.png, phi3.png, phi4.png с двумя кривыми
