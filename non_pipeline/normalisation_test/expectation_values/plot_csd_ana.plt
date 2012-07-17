#!/usr/bin/local/gnuplot -persist
set terminal pngcairo enhanced color

set output "csd_ana_at_analysis.png"
set title "Analytical Cross-spectral density at analysis"
set ylabel "CSD"
set xlabel "freq[Hz]"
#set logscale
#set yrange [ 1e-8 : 1e4 ]
#set xrange [ 1e-5 : 1e-1 ]
#set ytics 10
#set mytics 2
set key top left
#set grid ytics mytics xtics
plot "csd_ana.dat" using 1:2 with lines title 'ana. Re{csd} analysis' , \
"csd_ana.dat" using 1:3 with lines title 'ana. Im{csd} analysis' , \
"PSD.dat" index 1:1 using 1:6 with lines title "ana. Re{PAE}" , \
"PSD.dat" index 1:1 using 1:7 with lines title "ana. Im{PAE}"