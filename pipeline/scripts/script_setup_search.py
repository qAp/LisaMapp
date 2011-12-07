#!/usr/bin/env python
import os
import sys
import cPickle as cpkl



setup = {}

""" Where are the executables for analysis?"""
setup['execdir'] = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/pipeline-running-scripts/new_executables/'

""" Where are the time-series? """
setup['tsdir'] = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/runs/EccenctricInclined_eta0_0_xi0_0_sw_1_t0_0/noise/TDI/Michelson/G2/AET/simulation/by_me/stime_1.0_N_86400_for_days/data/'


""" How much do you want to scale the amplitude of the time-series by? """
setup['scale_ts'] = 1e21


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


""" doX : Maximum likelihood estimation of X, the dirty map """
X = {}
X['orfdir'] = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/myLISAmodule_stuff/EccentricInclined_eta0_0_xi0_0_sw_1_t0_0/tdiORF_SpHs/'
X['IJs'] = [ 'AA' ]# , 'AT' , 'EE' , 'ET' , 'TT' ]
X['days'] = range( 1 , 365 + 1 )
X['GWslopes'] = [ 0 ]
X['lmax'] = 20
X['flow'] , X['fhigh'] = 3e-4 , 3e-2


""" do_G : Maximum likelihood estimation of G, the Fisher Matrix """
G = {}
G['orfdir'] = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/myLISAmodule_stuff/EccentricInclined_eta0_0_xi0_0_sw_1_t0_0/tdiORF_SpHs/'
G['IJs'] = [ 'AA' ] #[ 'AA' , 'AE' , 'AT' , 'EE' , 'ET' , 'TT' ]
G['days'] = range( 1 , 365+1 )
G['GWslopes'] = [ 0 ]
G['lmax'] = 20
G['flow'] , G['fhigh'] = 3e-4 , 3e-2


""" do_sigma_avg : Calculate the sky-average standard deviation of Plm """
sigma_avg = {}
sigma_avg['IJs'] = [ 'AA' ]
sigma_avg['GWslopes'] = [ 0 ]
sigma_avg['lmax'] = 20
sigma_avg['regMethod'] , sigma_avg['regCutoff'] = 1 , 9.5e-4
sigma_avg['nlon'] , sigma_avg['nlat'] = 180 , 91



""" Post-processing """
postproc = {}
postproc['antennas'] = [ 'AA' ] #[ 'AA' , 'AE' , 'AT' , 'EE' , 'ET' , 'TT' ]
postproc['proj_X'] = True
postproc['proj_P'] , postproc['save_P'] = True , True
postproc['proj_Sg'] , postproc['save_Sg'] = True , True
postproc['proj_SNR'] , postproc['save_SNR'] = True , True
postproc['proj_stdP'] = True
postproc['lmax'] = 10
postproc['regMethod'] , postproc['regCutoff'] = 1 , 9.5e-4
postproc['nlat'] , postproc['nlon'] = 41 , 80






setup['psd'] = psd  
setup['csd'] = csd 
setup['avgpsd'] = avgpsd
setup['X'] = X
setup['G'] = G
setup['sigma_avg'] = sigma_avg
setup['postproc'] = postproc

#save this dictionary to disk
file = open( 'setup.pkl' , 'wb' )
cpkl.dump( setup , file , -1 )
file.close()










