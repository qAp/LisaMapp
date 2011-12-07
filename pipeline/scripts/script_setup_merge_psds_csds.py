#!/usr/bin/env python
import os
import sys
import cPickle as cpkl
import pprint 



#mergs = [ ( 'noise'/'signal' , psddir , normalisation ) , ... ]

mergs = [
    ( 'noise' ,
      '/gpfs1/JC0311443/workhere/stochasGW/Mapp/runs/EccenctricInclined_eta0_0_xi0_0_sw_1_t0_0/noise/TDI/Michelson/G2/AET/expecteds/psd_csd_for_days/csd/' ,
      1. ) ] + [
    ( 'signal' ,
      '/gpfs1/JC0311443/workhere/stochasGW/Mapp/runs/EccenctricInclined_eta0_0_xi0_0_sw_1_t0_0/signal/point_sources/lon_263_lat_-35_plus_P00_GWslope_0/TDI/Michelson/G2/AET/expecteds/psd_csd_for_days/csd/' ,
      1e-36 )  ] 


setup = {}
setup['execdir'] = ''
setup['days'] = [ 66 , 166 , 266 ]
setup['mergs'] = mergs
setup['number of batches'] = 1

file = open( 'setup_merge_psds_csds.pkl' , 'wb' ) ; cpkl.dump( setup , file , -1 ) ; file.close()


