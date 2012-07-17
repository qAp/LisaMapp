#!/usr/bin/local/gnuplot -persist
set terminal postscript enhanced color

set output "figures_sd_residuals/g00.eps"
set title "Residual"
set ylabel "(g00_cIJ - g00_ana) / g00_ana"
set xlabel "freq[Hz]"
plot "sd_residuals.dat" using 1:2 title "Re{g00}" , "sd_residuals.dat" using 1:3 title "Im{g00}"

set output "figures_sd_residuals/P12.eps"
set title "Residual"
set ylabel "(P12_cIJ - P12_ana) / P12_ana"
set xlabel "freq[Hz]"
plot "sd_residuals.dat" using 1:4 title "Re{P12}" , "sd_residuals.dat" using 1:5 title "Im{P12}"

set output "figures_sd_residuals/P11.eps"
set title "Residual"
set ylabel "(P11_cIJ - P11_ana) / P11_ana"
set xlabel "freq[Hz]"
plot "sd_residuals.dat" using 1:6 title "Re{P11}" , "sd_residuals.dat" using 1:7 title "Im{P11}" 

set output "figures_sd_residuals/P22.eps"
set title "Residual"
set ylabel "(P22_cIJ - P22_ana) / P22_ana"
set xlabel "freq[Hz]"
plot "sd_residuals.dat" using 1:8 title "Re{P22}" , "sd_residuals.dat" using 1:9 title "Im{P22}" 
