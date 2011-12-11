#!/usr/bin/env python
import os
import sys
import cPickle as cpkl



setup = {}

""" Where are the executables for analysis?"""
setup['execdir'] = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/labs/temp/analyse_uncorrelated_white_noises_same_white_signal/'

""" Where are the time-series? """
setup['tsdir'] = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/labs/sample_outputs/uncorrelated_white_noises_same_white_signal/stime_0.5_N_86400_for_days/data/'

""" Where do want to save the results? psd/, avgpsd/ and GW_slope_* etc. will be saved in this directory """
setup['workdir'] = 'analysis_seed_104_stime_0.5_hanning'


""" How much do you want to scale the amplitude of the time-series by? """
setup['scale_ts'] = 1.


""" plot_ts : Plot time-series"""
plot_ts = {}
plot_ts['days'] = [ 66 , 166 , 266 ]
plot_ts['scale_ts'] = 1e-20

""" dopsd : PSD estimation """
psd = {}
psd['days'] = range( 1 , 365 + 1 )
psd['segduration'] = 5760

""" docsd : CSD Estimation """
csd = {}
csd['days'] = range( 1 , 365 + 1 )
csd['segduration'] = 5760

""" do_avgpsd : Average PSD from adjacent days """
avgpsd = {}
avgpsd['days'] = range( 1 , 365 + 1 )


""" do_X : Maximum likelihood estimation of X, the dirty map """
X = {}
X['orfdir'] = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/labs/sample_glms/'
X['IJs'] = [ 'AE' ]# , 'AT' , 'ET' ] #[ 'AA' , 'AE' , 'AT' , 'EE' , 'ET' , 'TT' ]
X['days'] = range( 1 , 365 + 1 )
X['GWslopes'] = [ 0 ]
X['lmax'] = 0
X['window'] = 'hanning'
X['flow'] , X['fhigh'] = 1.74e-4 , 4.998e-1


""" do_G : Maximum likelihood estimation of G, the Fisher Matrix """
G = {}
G['orfdir'] = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/labs/sample_glms/'
G['IJs'] = [ 'AE' ]# , 'AT' , 'ET' ] #[ 'AA' , 'AE' , 'AT' , 'EE' , 'ET' , 'TT' ]
G['days'] = range( 1 , 365+1 )
G['GWslopes'] = [ 0 ]
G['lmax'] = 0
G['flow'] , G['fhigh'] = 1.74e-4 , 4.998e-1


""" do_S : The 'Psi' matrix in the bias of the clean map's covariance in the strong-signal limit (only available analytically, so use x_analyse_test.py) """
S = {}
S['csddir'] = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/labs/sample_outputs/uncorrelated_white_noises_same_white_signal/stime_0.5_N_86400_for_days/csd_ana/'
S['orfdir'] = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/labs/sample_glms/'
S['IJs'] = [ 'AE' ]# , 'AT' , 'ET' ]
S['days'] = range( 1 , 365+1 )
S['GWslopes'] = [ 0 ]
S['lmax'] = 0
S['flow'] , S['fhigh'] = 1.74e-4 , 4.998e-1


setup['plot_ts'] = plot_ts
setup['psd'] = psd  
setup['csd'] = csd 
setup['avgpsd'] = avgpsd
setup['X'] = X
setup['G'] = G
setup['S'] = S

#save this dictionary to disk
file = open( 'setup_analysis.pkl' , 'wb' )
cpkl.dump( setup , file , -1 )
file.close()










