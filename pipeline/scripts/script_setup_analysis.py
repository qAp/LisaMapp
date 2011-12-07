#!/usr/bin/env python
import os
import sys
import cPickle as cpkl



setup = {}

""" Where are the executables for analysis?"""
setup['execdir'] = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/pipeline-running-scripts/new_executables/'

""" Where are the time-series? """
setup['tsdir'] = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/runs/EccenctricInclined_eta0_0_xi0_0_sw_1_t0_0/signal/spherical_harmonics/Y_l0_m0_GWSlope_0/TDI/Michelson/G2/AET/simulation/stime_1.0_N_86400_for_days/data_r50/'

""" Where do want to save the results? psd/, avgpsd/ and GW_slope_* etc. will be saved in this directory """
setup['workdir'] = '/r50'


""" How much do you want to scale the amplitude of the time-series by? """
setup['scale_ts'] = 1e20


""" plot_ts : Plot time-series"""
plot_ts = {}
plot_ts['days'] = [ 66 , 166 , 266 ]
plot_ts['scale_ts'] = 1e-20

""" dopsd : PSD estimation """
psd = {}
psd['days'] = range( 1 , 365 + 1 )
psd['segduration'] = 2 * 60.**2

""" docsd : CSD Estimation """
csd = {}
csd['days'] = range( 1 , 365 + 1 )
csd['segduration'] = 2 * 60**2

""" do_avgpsd : Average PSD from adjacent days """
avgpsd = {}
avgpsd['days'] = range( 1 , 365 + 1 )


""" do_X : Maximum likelihood estimation of X, the dirty map """
X = {}
X['orfdir'] = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/myLISAmodule_stuff/EccentricInclined_eta0_0_xi0_0_sw_1_t0_0/tdiORF_SpHs/'
X['IJs'] = [ 'AE' , 'AT' , 'ET' ] #[ 'AA' , 'AE' , 'AT' , 'EE' , 'ET' , 'TT' ]
X['days'] = range( 1 , 365 + 1 )
X['GWslopes'] = [ 0 ]
X['lmax'] = 20
X['window'] = 'hanning'
X['flow'] , X['fhigh'] = 1.5e-4 , 4.998e-1


""" do_G : Maximum likelihood estimation of G, the Fisher Matrix """
G = {}
G['orfdir'] = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/myLISAmodule_stuff/EccentricInclined_eta0_0_xi0_0_sw_1_t0_0/tdiORF_SpHs/'
G['IJs'] = [ 'AE' , 'AT' , 'ET' ] #[ 'AA' , 'AE' , 'AT' , 'EE' , 'ET' , 'TT' ]
G['days'] = range( 1 , 365+1 )
G['GWslopes'] = [ 0 ]
G['lmax'] = 20
G['flow'] , G['fhigh'] = 1.5e-4 , 4.998e-1


""" do_S : The 'Psi' matrix in the bias of the clean map's covariance in the strong-signal limit (only available analytically, so use x_analyse_test.py) """
S = {}
S['orfdir'] = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/myLISAmodule_stuff/EccentricInclined_eta0_0_xi0_0_sw_1_t0_0/tdiORF_SpHs/'
S['IJs'] = [ 'AE' , 'AT' , 'ET' ]
S['days'] = [ 2 , 3 ]
S['GWslopes'] = [ 0 ]
S['lmax'] = 20
S['flow'] , S['fhigh'] = 1.5e-4 , 4.998e-1


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










