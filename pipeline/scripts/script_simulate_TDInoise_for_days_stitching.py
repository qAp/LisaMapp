#!/usr/bin/env python
import os
import sys
import glob
import numpy as np
import myUsefuls as mufls
import time
import subprocess

days = range( 1 , 5+1 )
stime = 0.5
seed = '1'
Tseg = 43200
tsdir = 'data'
Ntogroup = 2
tditype , tdigen , whichtdis = 'Michelson' , 'G2' , 'AET'
afterokpath = ''
jobidspath = ''




days_batches = mufls.divide_days_in_batches_reverse_pairup( days , Ntogroup )
Nb = len( days_batches )    

if afterokpath in ['','None',None,0] :
    afteroks = ['']
elif afterokpath not in glob.glob( afterokpath ) :
    raise Exception, 'The file you specified that contains IDs of jobs to be completed first cannot be found.'
else :
    file = open( afterokpath , 'r' ) ; lines = file.readlines() ; file.close()
    afteroks = [ line.rstrip() for line in lines ]

jobids = []
for b in range( Nb ) :
    batch = b + 1
    print "Preparing and submitting batch %d" % batch
    if days_batches[ b ] == [] :
        print 'No days in this batch.  Nothing to do...' ; sys.exit()
    else :
        commands = [ ('./x_simulate_TDInoise_for_days_stitching.py --tditype %s --tdigen %s --whichtdis %s -d%d --Tseg %f --stime %f --seed %s %s\n') %
                     ( tditype , tdigen , whichtdis , day , Tseg , stime , seed , tsdir )
                     for day in days_batches[ b ] ]
        submitname = 'x_simulate_TDInoise_for_days_stitching_b%03d' % batch
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
        p = subprocess.Popen( 'qsub -W depend=afterok:%s %s' %
                              ( ':'.join( afteroks ) , submitname + '.sub' ) ,
                              shell=True , stdout=subprocess.PIPE )
        jobids += [ p.communicate()[0].rstrip() ]
        
    time.sleep( 0.8 )

file = open( jobidspath , 'w' )
for k in range( len( jobids ) ) :
    print >> file , '%s' % jobids[k]
file.close()
