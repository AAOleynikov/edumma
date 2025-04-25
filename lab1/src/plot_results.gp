set terminal pngcairo size 800,600 enhanced font "Verdana,10"
set output "temperature.png"

# Получение статистики по температуре (столбец 3)
stats 'build/data.txt' using 3 nooutput

# Настройка палитры (9 цветов)
set palette defined ( \
    STATS_min                     "#0000FF",  \
    STATS_min + 1/8.0*(STATS_max-STATS_min) "#00B3FF", \
    STATS_min + 2/8.0*(STATS_max-STATS_min) "#00FFFF", \
    STATS_min + 3/8.0*(STATS_max-STATS_min) "#00FFB3", \
    STATS_min + 4/8.0*(STATS_max-STATS_min) "#00FF00", \
    STATS_min + 5/8.0*(STATS_max-STATS_min) "#B3FF00", \
    STATS_min + 6/8.0*(STATS_max-STATS_min) "#FFFF00", \
    STATS_min + 7/8.0*(STATS_max-STATS_min) "#FFA500", \
    STATS_max                     "#FF0000" \
)

# Настройка графика
set xlabel "X"
set ylabel "Y"
set size ratio -1
set grid
set colorbox

# Отрисовка прямоугольников (исправленный синтаксис)
plot 'build/data.txt' \
    using ($1 - $4/2.0):($2 - $5/2.0):4:5:3 \
    with boxxyerror \
    fill palette \
    fs solid 1.0 \
    notitle