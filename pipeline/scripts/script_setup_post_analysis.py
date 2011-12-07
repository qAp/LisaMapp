#!/usr/bin/env python
import os
import sys
import cPickle as cpkl



setup = {}


setup['execdir'] = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/pipeline-running-scripts/new_executables/'

setup['workdir'] = '/r50'

""" do_network_G : Calculate G for a network """
network_G = {}
network_G['networks'] = [ { 'GWslope':0 , 'IJs':['AE','AT'] } ,
                          { 'GWslope':0 , 'IJs':['AE','ET'] } ,
                          { 'GWslope':0 , 'IJs':['AT','ET'] } ,
                          { 'GWslope':0 , 'IJs':['AE','AT','ET'] } ]

""" do_network_X : Calculate X for a network """
network_X = {}
network_X['networks'] = [ { 'GWslope':0 , 'IJs':['AE','AT'] } ,
                          { 'GWslope':0 , 'IJs':['AE','ET'] } ,
                          { 'GWslope':0 , 'IJs':['AT','ET'] } ,
                          { 'GWslope':0 , 'IJs':['AE','AT','ET'] } ]

""" do_P : Calculate the clean map """
P = {}
P['IJs'] = [ 'AE' , 'AT' , 'ET' , 'AE_AT' , 'AE_ET' , 'AT_ET' , 'AE_AT_ET' ] 
P['GWslopes'] = [ 0 ]
P['regMethod'] , P['regCutoff'] = 0 , 1e-3 #1 , 9.5e-4
P['lmax'] , P['nlon'] , P['nlat'] = 0 , 180 , 91
P['mapnorm'] = 1.


""" do_sigma_avg : Calculate the sky-average standard deviation of Plm """
sigma_avg = {}
sigma_avg['IJs'] = [ 'AE' , 'AT' , 'ET' , 'AE_AT' , 'AE_ET' , 'AT_ET' , 'AE_AT_ET' ] 
sigma_avg['GWslopes'] = [ 0 ]
sigma_avg['regMethod'] , sigma_avg['regCutoff'] = 1 , 1e-3
sigma_avg['lmax'] , sigma_avg['nlon'] , sigma_avg['nlat'] = 20 , 180 , 91


""" do_sigmamap : Calculate the sigma map in the sky """
sigmamap = {}
sigmamap['IJs'] = [ 'AE'  , 'AT' , 'ET' , 'AE_AT' , 'AE_ET' , 'AT_ET' , 'AE_AT_ET' ] 
sigmamap['GWslopes'] = [ 0 ]
sigmamap['regMethod'] , sigmamap['regCutoff'] = 1 , 1e-3 #1 , 9.5e-4
sigmamap['lmax']  , sigmamap['nlon'] , sigmamap['nlat'] = 20 , 180 , 91
sigmamap['mapnorm'] = 1.


""" do_SNR_map : Calculate the SNR map """
SNRmap = {}
SNRmap['IJs'] = [ 'AE' ] #, 'AT' , 'ET' , 'AE_AT' , 'AE_ET' , 'AT_ET' , 'AE_AT_ET' ] 
SNRmap['GWslopes'] = [ 0 ]
SNRmap['regMethod'] , SNRmap['regCutoff'] = 1 , 1e-3 #1 , 9.5e-4
SNRmap['lmax'] , SNRmap['nlon'] , SNRmap['nlat'] = 15 , 180 , 91
SNRmap['mapnorm'] = 1.


""" do stdP : Calculate the standard deviation of Plm """
stdP = {}
stdP['IJs'] = [ 'AE' , 'AT' , 'ET' , 'AE_AT' , 'AE_ET' , 'AT_ET' , 'AE_AT_ET' ] 
stdP['GWslopes'] = [ 0 ]
stdP['regMethod'] , stdP['regCutoff'] = 0 , 1e-3 #1 , 9.5e-4
stdP['lmax'] = 0
stdP['Spath'] = 'default'


""" plot_X : Plot the dirty map X """
plot_X = {}
plot_X['IJs'] = [ 'AE'  , 'AT' , 'ET' , 'AE_AT' , 'AE_ET' , 'AT_ET' , 'AE_AT_ET' ] 
plot_X['GWslopes'] = [ 0 ]
plot_X['lmax'] , plot_X['nlon'] , plot_X['nlat'] = 15 , 180 , 91
plot_X['mapnorm'] = 1.


""" plot_P : Plot the clean map P """
plot_P = {}
plot_P['IJs'] = [ 'AE' , 'AT' , 'ET' , 'AE_AT' , 'AE_ET' , 'AT_ET' , 'AE_AT_ET' ] 
plot_P['GWslopes'] = [ 0 ]
plot_P['lmax'] , plot_P['nlon'] , plot_P['nlat'] = 0 , 180 , 91
plot_P['mapnorm'] = 1e-40


""" plot_sigma_avg : Plot the sky-average standard deviation of Plm """
plot_sigma_avg = {}
plot_sigma_avg['IJs'] = [ 'AE' ] # , 'AT' , 'ET' , 'AE_AT' , 'AE_ET' , 'AT_ET' , 'AE_AT_ET' ] 
plot_sigma_avg['GWslopes'] = [ 0 ]


""" plot_sigmamap : Plot sigmamap """
plot_sigmamap = {}
plot_sigmamap['IJs'] = [ 'AE' ] # , 'AT' , 'ET' , 'AE_AT' , 'AE_ET' , 'AT_ET' , 'AE_AT_ET' ] 
plot_sigmamap['GWslopes'] = [ 0 ]
plot_sigmamap['lmax'] , plot_sigmamap['nlon'] , plot_sigmamap['nlat'] = 15 , 180 , 91
plot_sigmamap['mapnorm'] = '1.'

""" plot_SNRmap : Plot SNRmap """
plot_SNRmap = {}
plot_SNRmap['IJs'] = [ 'AE' ] # , 'AT' , 'ET' , 'AE_AT' , 'AE_ET' , 'AT_ET' , 'AE_AT_ET' ] 
plot_SNRmap['GWslopes'] = [ 0 ]
plot_SNRmap['lmax'] , plot_SNRmap['nlon'] , plot_SNRmap['nlat'] = 15 , 180 , 91
plot_SNRmap['mapnorm'] = 1.

""" plot_stdP : Plot standard deviation of Plm """
plot_stdP = {}
plot_stdP['IJs'] = [ 'AE' , 'AT' , 'ET' , 'AE_AT' , 'AE_ET' , 'AT_ET' , 'AE_AT_ET' ] 
plot_stdP['GWslopes'] = [ 0 ]
plot_stdP['norm'] = 1.
plot_stdP['lmax'] = 15

""" plot_singular_values : Plot the singular values of the original Fisher matrix """
plot_singular_values = {}
plot_singular_values['IJs'] = [ 'AE' , 'AT' , 'ET' , 'AE_AT' , 'AE_ET' , 'AT_ET' , 'AE_AT_ET' ] 
plot_singular_values['GWslopes'] = [ 0 ]
plot_singular_values['lmax'] = 20




setup['P'] = P
setup['network_G'] = network_G
setup['network_X'] = network_X
setup['sigmamap'] = sigmamap
setup['SNRmap'] = SNRmap
setup['stdP'] = stdP
setup['sigma_avg'] = sigma_avg
setup['plot_X'] = plot_X
setup['plot_P'] = plot_P
setup['plot_sigmamap'] = plot_sigmamap
setup['plot_SNRmap'] = plot_SNRmap
setup['plot_stdP'] = plot_stdP
setup['plot_sigma_avg'] = plot_sigma_avg
setup['plot_singular_values'] = plot_singular_values



file = open( 'setup_post_analysis.pkl' , 'wb' ) ; cpkl.dump( setup , file , -1 ) ; file.close()
