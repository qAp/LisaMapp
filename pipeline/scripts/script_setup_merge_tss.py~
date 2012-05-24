#!/usr/bin/env python
import os
import sys
import cPickle as cpkl
import pprint 



#mergs = [ ( 'noise'/'signal' , tsddir , normalisation ) , ... ]

mergs = [
    ( 'noise' ,
      '/gpfs1/JC0311443/workhere/stochasGW/Mapp/runs/EccenctricInclined_eta0_0_xi0_0_sw_1_t0_0/noise/TDI/Michelson/G2/AET/simulation/by_me/stime_1.0_N_86400_for_days/data/',
      1. ) ] + [
    ( 'signal' ,
      '/gpfs1/JC0311443/workhere/stochasGW/Mapp/runs/EccenctricInclined_eta0_0_xi0_0_sw_1_t0_0/signal/point_sources/lon_263_lat_-35_plus_P00_GWslope_0/TDI/Michelson/G2/AET/simulation/stime_1.0_N_86400_for_days/data/',
      1e-18 ) ] 


setup = {}
setup['execdir'] = '~/workhere/stochasGW/Mapp/pipeline-running-scripts/new_executables/'
setup['days'] = range( 1 , 365+1 )
setup['mergs'] = mergs
setup['number of batches'] = 10

file = open( 'setup_merge_tss.pkl' , 'wb' ) ; cpkl.dump( setup , file , -1 ) ; file.close()


