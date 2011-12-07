#!/usr/bin/env python
import os
import sys
import glob as glob
import numpy as np
import cPickle as cpkl




execdir = '' #'/gpfs1/JC0311443/workhere/stochasGW/Mapp/pipeline-running-scripts/new_executables/'

os.system( 'cp %s .' % ( execdir + 'x_simulate_expected_signal_PSD.py' ) )

workdir = os.getcwd()  + '/'

file = open( 'setup.pkl' , 'rb' )
setup = cpkl.load( file )
file.close()

sourcename , Pdir , GWSpectralSlope = setup['sourcename'] , setup['Pdir'] , setup['GWSpectralSlope']
orfdir , tditype , tdigen = setup['orfdir'] , setup['tditype'] , setup['tdigen']
days = setup['days']
Nb = setup['number of batches']

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

    setupb = {}
    setupb['sourcename'] , setupb['Pdir'] , setupb['GWSpectralSlope'] = sourcename , Pdir , GWSpectralSlope
    setupb['orfdir'] , setupb['tditype'] , setupb['tdigen'] = orfdir , tditype , tdigen
    setupb['days'] = days_batches[ b ]

    setupname = 'setup_b%03d.pkl' % batch
    file = open( setupname , 'wb' )
    cpkl.dump( setupb , file , -1 )
    file.close()

    submitname = 'x_simulate_expected_signal_PSD_b%03d.sub' % batch
    file = open( submitname , 'w' )
    file.writelines( [ '#!/bin/bash\n' ,
                       '#PBS -N %s\n' % submitname ,
                       '#PBS -q compute\n' ,
                       '#PBS -j oe\n' ,
                       '#PBS -l nodes=1:ppn=1\n' ,
                       '#PBS -l walltime=2:00:00\n' ,
                       'cd $PBS_O_WORKDIR\n' ,
                       '\n' ,
                       './x_simulate_expected_signal_PSD.py %s' % setupname ] )
    file.close()

    os.system( 'qsub %s' % submitname )

