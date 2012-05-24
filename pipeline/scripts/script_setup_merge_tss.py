#!/usr/bin/env python
import os
import sys
import cPickle as cpkl
import pprint 



#mergs = [ ( 'noise'/'signal' , tsddir , normalisation ) , ... ]

mergs = [
    ( 'noise' ,
      '' , 
      1. ) ] + [
    ( 'signal' ,
      '' ,
      1. ) ] 


setup = {}
setup['execdir'] = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/Mapp_codes/pipeline/executables/'
setup['days'] = range( 1 , 365+1 )
setup['mergs'] = mergs
setup['number of batches'] = 61
setup['afterokpath'] = ''
setup['jobidspath'] = ''


file = open( 'setup_merge_tss.pkl' , 'wb' ) ; cpkl.dump( setup , file , -1 ) ; file.close()


