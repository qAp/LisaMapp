#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy.random as random
import optparse

usage = """
%prog SETUP_SIGNAL_SIMULATION.PKL\n
SETUP_SIGNAL_SIMULATION.PKL --- setup file for signal simulation ( created by script_setup_signal_simulation.py )
"""
parser = optparse.OptionParser( usage = usage )

( options , args ) = parser.parse_args()

if len( args ) < 1 :
    parser.error( 'You must specify a setup file for signal simulation!' )

setupname = args[0]
file = open( setupname , 'rb' ) ; setup = cpkl.load( file ) ; file.close()

execdir = setup['execdir']
tsdir = setup['tsdir']

if 'x_P_to_TS.py' not in glob.glob( 'x_P_to_TS.py' ) :
    os.system( 'cp %s .' % ( execdir + '/x_P_to_TS.py' ) )

workdir = os.getcwd() + '/'

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
    days = days_batches[ b ]
#    os.system( ( './x_P_to_TS.py ' + '-d%d '*len(days) + '--GWslope %d --tditype %s --tdigen %s --whichtdi %s --lmax %d --stime %f --f0 %f --df %f --Nf %d %s %s %s' )
#               % tuple( days + [ setup['GWslope'] , setup['tditype'] , setup['tdigen'] , setup['whichtdi'] , setup['lmax'] , setup['stime'] ,
#                                 setup['f0'] , setup['df'] , setup['Nf'] , setup['Ppath'] , setup['orfdir'] , tsdir ] ) )

    submitname = 'x_P_to_TS_b%03d.sub' % batch
    file = open( submitname , 'w' )
    file.writelines( [ '#!/bin/bash\n' ,
                       '#PBS -N %s\n' % submitname ,
                       '#PBS -q compute\n' ,
                       '#PBS -o x_P_to_TS_b%03d.out\n' % batch ,
                       '#PBS -j oe\n' ,
                       '#PBS -l nodes=1:ppn=1\n' ,
                       '#PBS -l walltime=10:00:00\n' ,
                       'cd $PBS_O_WORKDIR\n' ,
                       '\n' ,
                       ( './x_P_to_TS.py ' + '-d%d '*len(days) + '--seed %d --GWslope %d --tditype %s --tdigen %s --whichtdi %s --lmax %d --stime %f --f0 %f --df %f --Nf %d %s %s %s\n' )
               % tuple( days + [ setup['seed'] , setup['GWslope'] , setup['tditype'] , setup['tdigen'] , setup['whichtdi'] , setup['lmax'] , setup['stime'] ,
                                 setup['f0'] , setup['df'] , setup['Nf'] , setup['Ppath'] , setup['orfdir'] , tsdir ] ) ,
                       '\n' ,
                       'echo done' ] )
    file.close()
    file = open( 'x_P_to_TS_b%03d.out' % batch , 'w' ) ; file.write( 'dummy output' ) ; file.close()
    print 'Submitting batch %d' % batch
    os.system( 'qsub %s' % submitname )
    print 'done'

    

    
