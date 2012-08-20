#!/usr/bin/gnuplot -persist
set terminal pdfcairo enhanced color

#set logscale
set output "figures_seedsavg_daysavg/P11.pdf"
set title "PSD:11"
set ylabel "P11[s]"
set xlabel "freq[Hz]"
plot "sd_seedsavg_daysavg/avgsd.dat" index 0:0 using 1:2 with lines title "P11" , \
"expectation_values/PSD.dat" index 0:0 using 1:6 with lines title "ana.P11"

set output "figures_seedsavg_daysavg/P22.pdf"
set title "PSD:22"
set ylabel "P22[s]"
set xlabel "freq[Hz]"
plot "sd_seedsavg_daysavg/avgsd.dat" index 0:0 using 1:3 with lines title "P22" , \
"expectation_values/PSD.dat" index 4:4 using 1:6 with lines title "ana.P22"


set output "figures_seedsavg_daysavg/P33.pdf"
set title "PSD:33"
set ylabel "P33[s]"
set xlabel "freq[Hz]"
plot "sd_seedsavg_daysavg/avgsd.dat" index 0:0 using 1:4 with lines title "P33" , \
"expectation_values/PSD.dat" index 8:8 using 1:6 with lines title "ana.P33"

#unset logscale
set output "figures_seedsavg_daysavg/P12.pdf"
set title "CSD:12"
set ylabel "P12[s]"
set xlabel "freq[Hz]"
plot "sd_seedsavg_daysavg/avgsd.dat" index 1:1 using 1:2 with lines title "Re{P12}" , \
"sd_seedsavg_daysavg/avgsd.dat" index 1:1 using 1:3 with lines title "Im{P12}" , \
"expectation_values/PSD.dat" index 1:1 using 1:6 with lines title "ana. Re{P12}" , \
"expectation_values/PSD.dat" index 1:1 using 1:7 with lines title "ana. Im{P12}"

set output "figures_seedsavg_daysavg/P13.pdf"
set title "CSD:13"
set ylabel "P13[s]"
set xlabel "freq[Hz]"
plot "sd_seedsavg_daysavg/avgsd.dat" index 1:1 using 1:4 with lines title "Re{P13}" , \
"sd_seedsavg_daysavg/avgsd.dat" index 1:1 using 1:5 with lines title "Im{P13}" , \
"expectation_values/PSD.dat" index 2:2 using 1:6 with lines title "ana. Re{P13}" , \
"expectation_values/PSD.dat" index 2:2 using 1:7 with lines title "ana. Im{P13}"

set output "figures_seedsavg_daysavg/P23.pdf"
set title "CSD:23"
set ylabel "P23[s]"
set xlabel "freq[Hz]"
plot "sd_seedsavg_daysavg/avgsd.dat" index 1:1 using 1:6 with lines title "Re{P23}" , \
"sd_seedsavg_daysavg/avgsd.dat" index 1:1 using 1:7 with lines title "Im{P23}" , \
"expectation_values/PSD.dat" index 5:5 using 1:6 with lines title "ana. Re{P23}" , \
"expectation_values/PSD.dat" index 5:5 using 1:7 with lines title "ana. Im{P23}"
