#!/usr/bin/local/gnuplot -persist

set terminal pngcairo enhanced color

set output "figures_PSDs/P11.png"
set title "Power spectral densities"
set ylabel "PSD"
set xlabel "freq[Hz]"
set key left top
plot "PSD.dat" index 0:0 using 1:2 with lines title "RePhh" , \
"PSD.dat" index 0:0 using 1:3 with lines title "ImPhh" , \
"PSD.dat" index 0:0 using 1:4 with lines title "RePnn" , \
"PSD.dat" index 0:0 using 1:5 with lines title "ImPnn" , \
"PSD.dat" index 0:0 using 1:6 with lines title "RePoo" , \
"PSD.dat" index 0:0 using 1:7 with lines title "ImPoo"

set output "figures_PSDs/P12.png"
set title "Power spectral densities"
set ylabel "PSD"
set xlabel "freq[Hz]"
set key left top
plot "PSD.dat" index 1:1 using 1:2 with lines title "RePhh" , \
"PSD.dat" index 1:1 using 1:3 with lines title "ImPhh" , \
"PSD.dat" index 1:1 using 1:4 with lines title "RePnn" , \
"PSD.dat" index 1:1 using 1:5 with lines title "ImPnn" , \
"PSD.dat" index 1:1 using 1:6 with lines title "RePoo" , \
"PSD.dat" index 1:1 using 1:7 with lines title "ImPoo"

set output "figures_PSDs/P13.png"
set title "Power spectral densities"
set ylabel "PSD"
set xlabel "freq[Hz]"
set key left top
plot "PSD.dat" index 2:2 using 1:2 with lines title "RePhh" , \
"PSD.dat" index 2:2 using 1:3 with lines title "ImPhh" , \
"PSD.dat" index 2:2 using 1:4 with lines title "RePnn" , \
"PSD.dat" index 2:2 using 1:5 with lines title "ImPnn" , \
"PSD.dat" index 2:2 using 1:6 with lines title "RePoo" , \
"PSD.dat" index 2:2 using 1:7 with lines title "ImPoo"

set output "figures_PSDs/P21.png"
set title "Power spectral densities"
set ylabel "PSD"
set xlabel "freq[Hz]"
set key left top
plot "PSD.dat" index 3:3 using 1:2 with lines title "RePhh" , \
"PSD.dat" index 3:3 using 1:3 with lines title "ImPhh" , \
"PSD.dat" index 3:3 using 1:4 with lines title "RePnn" , \
"PSD.dat" index 3:3 using 1:5 with lines title "ImPnn" , \
"PSD.dat" index 3:3 using 1:6 with lines title "RePoo" , \
"PSD.dat" index 3:3 using 1:7 with lines title "ImPoo"

set output "figures_PSDs/P22.png"
set title "Power spectral densities"
set ylabel "PSD"
set xlabel "freq[Hz]"
set key left top
plot "PSD.dat" index 4:4 using 1:2 with lines title "RePhh" , \
"PSD.dat" index 4:4 using 1:3 with lines title "ImPhh" , \
"PSD.dat" index 4:4 using 1:4 with lines title "RePnn" , \
"PSD.dat" index 4:4 using 1:5 with lines title "ImPnn" , \
"PSD.dat" index 4:4 using 1:6 with lines title "RePoo" , \
"PSD.dat" index 4:4 using 1:7 with lines title "ImPoo"

set output "figures_PSDs/P23.png"
set title "Power spectral densities"
set ylabel "PSD"
set xlabel "freq[Hz]"
set key left top
plot "PSD.dat" index 5:5 using 1:2 with lines title "RePhh" , \
"PSD.dat" index 5:5 using 1:3 with lines title "ImPhh" , \
"PSD.dat" index 5:5 using 1:4 with lines title "RePnn" , \
"PSD.dat" index 5:5 using 1:5 with lines title "ImPnn" , \
"PSD.dat" index 5:5 using 1:6 with lines title "RePoo" , \
"PSD.dat" index 5:5 using 1:7 with lines title "ImPoo"

set output "figures_PSDs/P31.png"
set title "Power spectral densities"
set ylabel "PSD"
set xlabel "freq[Hz]"
set key left top
plot "PSD.dat" index 6:6 using 1:2 with lines title "RePhh" , \
"PSD.dat" index 6:6 using 1:3 with lines title "ImPhh" , \
"PSD.dat" index 6:6 using 1:4 with lines title "RePnn" , \
"PSD.dat" index 6:6 using 1:5 with lines title "ImPnn" , \
"PSD.dat" index 6:6 using 1:6 with lines title "RePoo" , \
"PSD.dat" index 6:6 using 1:7 with lines title "ImPoo"

set output "figures_PSDs/P32.png"
set title "Power spectral densities"
set ylabel "PSD"
set xlabel "freq[Hz]"
set key left top
plot "PSD.dat" index 7:7 using 1:2 with lines title "RePhh" , \
"PSD.dat" index 7:7 using 1:3 with lines title "ImPhh" , \
"PSD.dat" index 7:7 using 1:4 with lines title "RePnn" , \
"PSD.dat" index 7:7 using 1:5 with lines title "ImPnn" , \
"PSD.dat" index 7:7 using 1:6 with lines title "RePoo" , \
"PSD.dat" index 7:7 using 1:7 with lines title "ImPoo"

set output "figures_PSDs/P33.png"
set title "Power spectral densities"
set ylabel "PSD"
set xlabel "freq[Hz]"
set key left top
plot "PSD.dat" index 8:8 using 1:2 with lines title "RePhh" , \
"PSD.dat" index 8:8 using 1:3 with lines title "ImPhh" , \
"PSD.dat" index 8:8 using 1:4 with lines title "RePnn" , \
"PSD.dat" index 8:8 using 1:5 with lines title "ImPnn" , \
"PSD.dat" index 8:8 using 1:6 with lines title "RePoo" , \
"PSD.dat" index 8:8 using 1:7 with lines title "ImPoo"




