#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import myLISAmodule as mlisar
import AnisotropySearch as AS


days = range( 1 , 365+1 ); #days.remove( 1 ) ; days.remove( 180 ) ; days.remove( 270 )


dayinsecs = 86400

tditype , tdigen , whichtdis = 'Michelson' , 'G2' , 'AET'
stime , duration = 1. , dayinsecs


for day in days :
    print 'Generating the noise for day %d' % day
    
    inittime = dayinsecs * ( day - 1 )

    tspath = 'data_r2/d%03d.pkl' % day
    setup = { 'tditype':tditype , 'tdigen':tdigen , 'whichtdis':whichtdis ,
              'stime':stime , 'duration':duration , 'inittime':inittime , 'seed':None ,
              'tspath':tspath }

    setupname = 'setup_d%03d.pkl' % day
    file = open( setupname , 'wb' ) ; cpkl.dump( setup , file , -1 ) ; file.close()

    submitname = 'd%03d.sub' % day
    file = open( submitname , 'w' )
    file.writelines( [ '#!/bin/bash\n' ,
                       '#PBS -N %s\n' % submitname ,
                       '#PBS -q compute\n' ,
                       '#PBS -j oe\n' ,
                       '#PBS -l nodes=1:ppn=1\n' ,
                       '#PBS -l walltime=1:00:00\n' ,
                       'cd $PBS_O_WORKDIR\n' ,
                       '\n' ,
                       './x_simulate_noise.py %s' % setupname ] ) ; file.close()

    os.system( 'qsub %s' % submitname )
