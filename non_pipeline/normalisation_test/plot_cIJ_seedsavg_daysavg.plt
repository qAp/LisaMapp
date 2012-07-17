#!/usr/bin/gnuplot -persist
set terminal postscript enhanced color

set output "figures_seedsavg_daysavg_cIJ/P11.eps"
set title "PSD:11"
set ylabel "P11[s]"
set xlabel "freq[Hz]"
plot "sd_seedsavg_daysavg_cIJ/avgsd.dat" index 0:0 using 1:6 with dots title "Re{P11}" , \
"sd_seedsavg_daysavg_cIJ/avgsd.dat" index 0:0 using 1:7 with dots title "Im{P11}" , \
"expectation_values/cIJ_ana/anasd.dat" using 1:6 with dots title "ana.Re{P11}" , \
"expectation_values/cIJ_ana/anasd.dat" using 1:7 with dots title "ana.Im{P11}"
#"expectation_values/PSD.dat" index 0:0 using 1:6 with lines title "ana.P11"

set output "figures_seedsavg_daysavg_cIJ/P22.eps"
set title "PSD:22"
set ylabel "P22[s]"
set xlabel "freq[Hz]"
plot "sd_seedsavg_daysavg_cIJ/avgsd.dat" index 0:0 using 1:8 with dots title "Re{P22}" , \
"sd_seedsavg_daysavg_cIJ/avgsd.dat" index 0:0 using 1:9 with dots title "Im{P22}" , \
"expectation_values/cIJ_ana/anasd.dat" using 1:8 with dots title "ana.Re{P22}" , \
"expectation_values/cIJ_ana/anasd.dat" using 1:9 with dots title "ana.Im{P22}" 
#"expectation_values/PSD.dat" index 4:4 using 1:6 with lines title "ana.P22"

set output "figures_seedsavg_daysavg_cIJ/P12.eps"
set title "CSD:12"
set ylabel "P12[s]"
set xlabel "freq[Hz]"
plot "sd_seedsavg_daysavg_cIJ/avgsd.dat" index 0:0 using 1:4 with dots title "Re{P12}" , \
"sd_seedsavg_daysavg_cIJ/avgsd.dat" index 0:0 using 1:5 with dots title "Im{P12}" , \
"expectation_values/cIJ_ana/anasd.dat" using 1:4 with dots title "ana.Re{P12}" , \
"expectation_values/cIJ_ana/anasd.dat" using 1:5 with dots title "ana.Im{P12}"
#"expectation_values/PSD.dat" index 1:1 using 1:6 with lines title "ana. Re{P12}" , \
#"expectation_values/PSD.dat" index 1:1 using 1:7 with lines title "ana. Im{P12}"

set output "figures_seedsavg_daysavg_cIJ/g12.eps"
set title "ORF_SpH: g12"
set ylabel "g12[s]"
set xlabel "freq[Hz]"
plot "sd_seedsavg_daysavg_cIJ/avgsd.dat" index 0:0 using 1:2 with dots title "Re{g12}" , \
"sd_seedsavg_daysavg_cIJ/avgsd.dat" index 0:0 using 1:3 with dots title "Im{g12}" , \
"expectation_values/cIJ_ana/anasd.dat" using 1:2 with dots title "ana.Re{g12}" , \
"expectation_values/cIJ_ana/anasd.dat" using 1:3 with dots title "ana.Im{g12}"

#set output "figures_seedsavg_daysavg_cIJ/P13.eps"
#set title "CSD:13"
#set ylabel "P13[s]"
#set xlabel "freq[Hz]"
#plot "sd_seedsavg_daysavg_cIJ/avgsd.dat" index 1:1 using 1:4 with lines title "Re{P13}" , \
#"sd_seedsavg_daysavg_cIJ/avgsd.dat" index 1:1 using 1:5 with lines title "Im{P13}" , \
#"expectation_values/PSD.dat" index 2:2 using 1:6 with lines title "ana. Re{P13}" , \
#"expectation_values/PSD.dat" index 2:2 using 1:7 with lines title "ana. Im{P13}"
#
#set output "figures_seedsavg_daysavg_cIJ/P23.eps"
#set title "CSD:23"
#set ylabel "P23[s]"
#set xlabel "freq[Hz]"
#plot "sd_seedsavg_daysavg_cIJ/avgsd.dat" index 1:1 using 1:6 with lines title "Re{P23}" , \
#"sd_seedsavg_daysavg_cIJ/avgsd.dat" index 1:1 using 1:7 with lines title "Im{P23}" , \
#"expectation_values/PSD.dat" index 5:5 using 1:6 with lines title "ana. Re{P23}" , \
#"expectation_values/PSD.dat" index 5:5 using 1:7 with lines title "ana. Im{P23}"
