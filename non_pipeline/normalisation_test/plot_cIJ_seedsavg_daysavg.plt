#!/usr/bin/gnuplot -persist
set terminal pdfcairo enhanced color

set arrow 1 from 0.0005, graph 0 to 0.0005, graph 1 nohead
set arrow 2 from 0.0620, graph 0 to 0.0620, graph 1 nohead

set output "figures_cIJ_seedsavg_daysavg/P11.pdf"
set title "PSD:11"
set ylabel "P11[s]"
set xlabel "freq[Hz]"
plot "cIJ_seedsavg_daysavg/avgsd.dat" index 0:0 using 1:6 with dots title "Re{P11}" , \
"cIJ_seedsavg_daysavg/avgsd.dat" index 0:0 using 1:7 with dots title "Im{P11}" , \
"expectation_values/cIJ_ana/anasd.dat" using 1:6 with dots title "ana.Re{P11}" , \
"expectation_values/cIJ_ana/anasd.dat" using 1:7 with dots title "ana.Im{P11}"
#"expectation_values/PSD.dat" index 0:0 using 1:6 with lines title "ana.P11"

set output "figures_cIJ_seedsavg_daysavg/P22.pdf"
set title "PSD:22"
set ylabel "P22[s]"
set xlabel "freq[Hz]"
plot "cIJ_seedsavg_daysavg/avgsd.dat" index 0:0 using 1:8 with dots title "Re{P22}" , \
"cIJ_seedsavg_daysavg/avgsd.dat" index 0:0 using 1:9 with dots title "Im{P22}" , \
"expectation_values/cIJ_ana/anasd.dat" using 1:8 with dots title "ana.Re{P22}" , \
"expectation_values/cIJ_ana/anasd.dat" using 1:9 with dots title "ana.Im{P22}" 
#"expectation_values/PSD.dat" index 4:4 using 1:6 with lines title "ana.P22"

set output "figures_cIJ_seedsavg_daysavg/P12.pdf"
set title "CSD:12"
set ylabel "P12[s]"
set xlabel "freq[Hz]"
plot "cIJ_seedsavg_daysavg/avgsd.dat" index 0:0 using 1:4 with dots title "Re{P12}" , \
"cIJ_seedsavg_daysavg/avgsd.dat" index 0:0 using 1:5 with dots title "Im{P12}" , \
"expectation_values/cIJ_ana/anasd.dat" using 1:4 with dots title "ana.Re{P12}" , \
"expectation_values/cIJ_ana/anasd.dat" using 1:5 with dots title "ana.Im{P12}"
#"expectation_values/PSD.dat" index 1:1 using 1:6 with lines title "ana. Re{P12}" , \
#"expectation_values/PSD.dat" index 1:1 using 1:7 with lines title "ana. Im{P12}"

set output "figures_cIJ_seedsavg_daysavg/g12.pdf"
set title "ORF_SpH: g12"
set ylabel "g12[s]"
set xlabel "freq[Hz]"
plot "cIJ_seedsavg_daysavg/avgsd.dat" index 0:0 using 1:2 with dots title "Re{g12}" , \
"cIJ_seedsavg_daysavg/avgsd.dat" index 0:0 using 1:3 with dots title "Im{g12}" , \
"expectation_values/cIJ_ana/anasd.dat" using 1:2 with dots title "ana.Re{g12}" , \
"expectation_values/cIJ_ana/anasd.dat" using 1:3 with dots title "ana.Im{g12}"

#set output "figures_cIJ_seedsavg_daysavg/P13.pdf"
#set title "CSD:13"
#set ylabel "P13[s]"
#set xlabel "freq[Hz]"
#plot "cIJ_seedsavg_daysavg/avgsd.dat" index 1:1 using 1:4 with lines title "Re{P13}" , \
#"cIJ_seedsavg_daysavg/avgsd.dat" index 1:1 using 1:5 with lines title "Im{P13}" , \
#"expectation_values/PSD.dat" index 2:2 using 1:6 with lines title "ana. Re{P13}" , \
#"expectation_values/PSD.dat" index 2:2 using 1:7 with lines title "ana. Im{P13}"
#
#set output "figures_cIJ_seedsavg_daysavg/P23.pdf"
#set title "CSD:23"
#set ylabel "P23[s]"
#set xlabel "freq[Hz]"
#plot "cIJ_seedsavg_daysavg/avgsd.dat" index 1:1 using 1:6 with lines title "Re{P23}" , \
#"cIJ_seedsavg_daysavg/avgsd.dat" index 1:1 using 1:7 with lines title "Im{P23}" , \
#"expectation_values/PSD.dat" index 5:5 using 1:6 with lines title "ana. Re{P23}" , \
#"expectation_values/PSD.dat" index 5:5 using 1:7 with lines title "ana. Im{P23}"
