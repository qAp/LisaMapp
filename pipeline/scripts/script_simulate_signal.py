#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy.random as random


#Where is the executable for this?
execdir = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/pipeline-running-scripts/new_executables/'

if 'x_P_to_TS.py' not in glob.glob( 'x_P_to_TS.py' ) :
    os.system( 'cp %s .' % ( execdir + 'x_P_to_TS.py' ) )

workdir = os.getcwd() + '/'

file = open( 'setup.pkl' , 'rb' )
setup = cpkl.load( file )
file.close()

sourcename = setup['sourcename']
Pdir = setup['Pdir']
GWSpectralSlope = setup['GWSpectralSlope']

orfdir = setup['orfdir'] 
tditype , tdigen = setup['tditype'] , setup['tdigen']
lmax = setup['lmax']
f0 , df , Nf = setup['f0'] , setup['df'] , setup['Nf']

days = setup['days']
Nb = setup['number of batches']
seed = setup['seed']


random.seed( seed )

Ndays = len( days )
if Ndays % Nb == 0 :
    Ndb = Ndays / Nb
    days_batches = [ days[ b*Ndb : (b+1)*Ndb  ]  for b in range( Nb ) ]
elif Ndays % Nb > 0 :
    Ndb = Ndays / ( Nb-1 )
    Ndbl = Ndays % ( Nb-1 )
    days_batches = [ days[ b*Ndb : (b+1)*Ndb ] for b in range( Nb-1 ) ] + [ days[ (Nb-1)*Ndb : ] ]
else :
    raise Exception , "Both the number of days and number of batches have to be postivie integer!"


for b in range( Nb ) :

    batch = b + 1
    print "~~~~~~~ Processing batch %d ~~~~~~~" % batch
    
    setup = {}
    setup['sourcename'] = sourcename
    setup['Pdir'] = Pdir
    setup['GWSpectralSlope'] = GWSpectralSlope

    setup['orfdir'] = orfdir
    setup['tditype'] , setup['tdigen'] = tditype , tdigen
    setup['lmax'] = lmax
    setup['f0'] , setup['df'] , setup['Nf'] = f0 , df , Nf

    setup['days'] = days_batches[ b ]

    setupname = 'setup_b%03d.pkl' % batch
    file = open( setupname , 'wb' )
    cpkl.dump( setup , file , -1 )
    file.close()

    print "batch %d is set up, writing submit file..." % batch

    submitname = 'x_P_to_TS_b%03d.sub' % batch
    file = open( submitname , 'w' )
    lines = [ '#!/bin/bash\n' ,
              '#PBS -N %s\n' % submitname ,
              '#PBS -q compute\n' ,
              '#PBS -j oe\n' ,
              '#PBS -l nodes=1:ppn=1\n' ,
              '#PBS -l walltime=5:00:00\n' ,
              'cd $PBS_O_WORKDIR\n' ,
              '\n' ,
              './x_P_to_TS.py %s' % setupname ]
    file.writelines( lines )
    file.close()

    os.system( 'qsub %s' % submitname )
    

    

    
