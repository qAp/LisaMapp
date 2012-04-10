#!/usr/bin/env python
import os
import sys
import glob
import numpy as np
import myUsefuls as mufls
import time


days = range( 1 , 5+1 )
stime = 0.5
seed = 100
Tseg = 43200
GWSpectralSlope = 0
lmax = 0
compute_ORF_SpHs = True
Ppath = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/skymaps/library/sphericalhar\
monics/Y_l0_m0_x1_over_90_lmax_0/Y_l0_m0.pkl'
orfdir = 'temp_orfs'
tsdir = 'data'
Ntogroup = 2



days_batches = mufls.divide_days_in_batches_reverse_pairup( days , Ntogroup )
Nb = len( days_batches )    

for b in range( Nb ) :
    batch = b + 1
    print "Preparing and submitting batch %d" % batch
    if days_batches[ b ] == [] :
        print 'No days in this batch.  Nothing to do...' ; sys.exit()
    else :
        commands = [ ('./x_simulate_AETnoise_from_arbitrary_SpH_for_days_stitching.py \
        -d%d --Tseg %f --stime %f --GWSpectralSlope %d --lmax %d --seed %s \
        --compute_ORF_SpHs %d %s %s %s\n') %
                     ( day , Tseg , stime , GWSpectralSlope , lmax , seed ,\
                       compute_ORF_SpHs , orfdir , Ppath , tsdir )
                     for day in days_batches[ b ] ]
        submitname = 'x_simulate_AETnoiseX_b%03d' % batch
        file = open( submitname + '.sub' , 'w' )
        file.writelines( [ '#!/bin/bash\n' ,
                           '#PBS -N %s\n' % ( submitname + '.sub' ) ,
                           '#PBS -o %s\n' % ( submitname + '.out' ) ,
                           '#PBS -j oe\n' ,
                           '#PBS -q compute\n' ,
                           '#PBS -l nodes=1:ppn=1\n' ,
                           '#PBS -l walltime=10:00:00\n' , 
                           'cd $PBS_O_WORKDIR\n' , '\n' ] +
                         commands + [ 'echo done' ] ) ; file.close()
        file = open( submitname + '.out' , 'w' ) ; file.write( 'dummy output' ) ; file.close()
        os.system( 'qsub %s' % ( submitname + '.sub' ) ) ; print 'done'
    time.sleep( 0.5 )

