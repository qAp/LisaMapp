#!/usr/local/bin/gnuplot -persist
set terminal pngcairo color enhanced


set output "stats_mean.png"
set title "Sample mean of time-series"
set ylabel "{/Symbol m}"
set xlabel "time [day]"
plot "stats.dat" using 1:2 every :::::2 with lines title "{/Symbol m}_{1}" , \
"stats.dat" using 1:3 every :::::3 with lines title "{/Symbol m}_{2}" ,\
0 


set output "stats_variance.png"
set title "Sample variance of time-series"
set ylabel "{/Symbol s}"
set xlabel "time [day]"
plot "stats.dat" using 1:4 every :::::4 with lines title "{/Symbol s}_{11}" ,\
"stats.dat" using 1:7 every :::::7 with lines title "{/Symbol s}_{22}" , \
0.00987723338485 , \
0.0102396679404


set output "stats_covariance.png"
set title "Sample covariance of time-series"
set ylabel "{/Symbol s}"
set xlabel "time [day]"
plot "stats.dat" using 1:5 every :::::5 with lines title "{/Symbol s}_{12}" ,\
"stats.dat" using 1:6 every :::::6 with lines title "{/Symbol s}_{21}" , \
-1.88203602997e-05
