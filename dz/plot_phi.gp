# plot_phi_selected.gp
# Ожидается: datafile с колонками "# idx t phi1 phi3 phi4"

datafile = "data.dat"
img_w = 1200
img_h = 800

set datafile separator whitespace
set format x "%g"
set format y "%g"
set xlabel "t (s)"
set ylabel "phi (units)"
set grid
set style line 1 lw 2 pointsize 0.5

# Массив колонок: phi1 -> col 3, phi3 -> col 4, phi4 -> col 5
# Соответствие имён для файлов и заголовков
names = "phi1 phi3 phi4"
cols  = "3 4 5"

do for [i=1:3] {
    # получаем номер столбца и имя
    col = int(word(cols, i))
    nm  = word(names, i)
    outfile = sprintf("%s.png", nm)
    titletext = sprintf("%s (column %d)", nm, col)

    set terminal pngcairo size img_w,img_h enhanced font "Arial,12"
    set output outfile

    set title titletext
    plot datafile using 2:col with linespoints ls 1 title titletext

    set output
}

# Готово — в каталоге появятся файлы: phi1.png, phi3.png, phi4.png
