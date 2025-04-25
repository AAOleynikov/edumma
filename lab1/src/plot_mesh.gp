set terminal pngcairo size 4000,4000 enhanced font "Arial,10"
set output 'mesh.png'

# Чтение CSV с комментариями и разделителем
set datafile separator ","
set datafile commentschars "#"

unset key           # Убираем легенду
unset colorbox      # Убираем шкалу цветов

set xlabel 'x'
set ylabel 'y'
set grid
set size ratio 1

# Определение палитры для цветов точек по bc_type
set palette defined ( 0 "blue", 1 "green", 2 "red", 3 "purple" )

# Отрисовка:
#   точки увеличены (ps 1.5), подписи i,j и bc_val мелким шрифтом
plot \
    'build/mesh_nodes.csv' using 3:4:5 with points pt 7 ps 6 palette notitle, \
    ''                  using 3:4:(sprintf("%d,%d val=%.2f", column(1), column(2), column(6))) \
                       with labels font "Arial,8" offset char 0.5,0.5 center notitle