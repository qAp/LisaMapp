#!/usr/bin/env python
import os
import sys
import glob as glob
import numpy as np
import cPickle as cpkl
import optparse

usage = """
%prog SETUP_SIMULATE_EXPECTED_SIGNAL_PSD_OR_CSD.PKL\n
SETUP_SIMULATE_EXPECTED_SIGNAL_PSD_OR_CSD.PKL --- setup file ( generated by script_setup_simulate_expected_signal_PSD_or_CSD.py )
"""

parser = optparse.OptionParser( usage = usage )

parser.add_option( '--do_expected_PSDs' , action='store_true' , help='Get theoretical PSDs' )

parser.add_option( '--do_expected_CSDs' , action='store_true' , help='Get theoretical CSDs' )

( options , args ) = parser.parse_args()

if len( args ) < 1 :
    parser.error( 'You must specify at least SETUP_SIMULATE_EXPECTED_SIGNAL_PSD_OR_CSD.PKL!' )
    
file = open( args[ 0 ] , 'rb' ) ; setup = cpkl.load( file ) ; file.close()

workdir = os.getcwd() + '/'

Nb = setup['number of batches']

Ndays = len( setup['days'] )
if Ndays % Nb == 0 :
    Ndb = Ndays / Nb
    days_batches = [ setup['days'][ b*Ndb : (b+1)*Ndb  ]  for b in range( Nb ) ]
elif Ndays % Nb > 0 :
    Ndb = Ndays / ( Nb-1 )
    Ndbl = Ndays % ( Nb-1 )
    days_batches = [ setup['days'][ b*Ndb : (b+1)*Ndb ] for b in range( Nb-1 ) ] + [ setup['days'][ (Nb-1)*Ndb : ] ]
else :
    raise Exception , "Both the number of days and number of batches have to be postivie integer!"


if options.do_expected_PSDs :
    print 'Getting the theoretical PSDs...'
    if 'x_simulate_expected_signal_PSD.py' not in glob.glob( 'x_simulate_expected_signal_PSD.py' ) :
        os.system( 'cp %s .' % ( setup['execdir'] + 'x_simulate_expected_signal_PSD.py' ) )

    psddir = workdir + 'psd/'

    for b in range( Nb ) :
        batch = b + 1 ; print "Processing batch %d " % batch
        days = days_batches[ b ]
        print days

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
                           ( './x_simulate_expected_signal_PSD.py ' + '-d%d '*len(days) + '--GWslope %d --tditype %s --tdigen %s --whichtdi %s --lmax %d --f0 %f --df %f --Nf %d %s %s %s ' ) %
                           tuple( days + [ setup['GWslope'] , setup['tditype'] , setup['tdigen'] , setup['whichtdi'] ,
                                           setup['lmax'] , setup['f0'] , setup['df'] , setup['Nf'] ,
                                           setup['Ppath'] , setup['orfdir'] , psddir ] ) ] )
        file.close()
        os.system( 'qsub %s' % submitname )
        print 'done'
#        os.system( ( './x_simulate_expected_signal_PSD.py ' + '-d%d '*len(days) + '--GWslope %d --tditype %s --tdigen %s --whichtdi %s --lmax %d --f0 %f --df %f --Nf %d %s %s %s ' ) %
#                   tuple( days + [ setup['GWslope'] , setup['tditype'] , setup['tdigen'] , setup['whichtdi'] ,
#                                            setup['lmax'] , setup['f0'] , setup['df'] , setup['Nf'] ,
#                                            setup['Ppath'] , setup['orfdir'] , psddir ] ) )
    print 'done'




if options.do_expected_CSDs :
    print 'Getting the theoretical CSDs...'
    if 'x_simulate_expected_signal_CSD.py' not in glob.glob( 'x_simulate_expected_signal_CSD.py' ) :
        os.system( 'cp %s .' % ( setup['execdir'] + 'x_simulate_expected_signal_CSD.py' ) )

    csddir = workdir + 'csd/'

    for b in range( Nb ) :
        batch = b + 1 ; print "Processing batch %d " % batch
        days = days_batches[ b ]
        print days

        submitname = 'x_simulate_expected_signal_CSD_b%03d.sub' % batch
        file = open( submitname , 'w' )
        file.writelines( [ '#!/bin/bash\n' ,
                           '#PBS -N %s\n' % submitname ,
                           '#PBS -q compute\n' ,
                           '#PBS -j oe\n' ,
                           '#PBS -l nodes=1:ppn=1\n' ,
                           '#PBS -l walltime=2:00:00\n' ,
                           'cd $PBS_O_WORKDIR\n' ,
                           '\n' ,
                           ( './x_simulate_expected_signal_CSD.py ' + '-d%d '*len(days) + '--GWslope %d --tditype %s --tdigen %s --whichtdi %s --lmax %d --f0 %f --df %f --Nf %d %s %s %s ' ) %
                           tuple( days + [ setup['GWslope'] , setup['tditype'] , setup['tdigen'] , setup['whichtdi'] ,
                                           setup['lmax'] , setup['f0'] , setup['df'] , setup['Nf'] ,
                                           setup['Ppath'] , setup['orfdir'] , csddir ] ) ] )
        file.close()
        os.system( 'qsub %s' % submitname )        
        print 'done'
#        os.system( ( './x_simulate_expected_signal_CSD.py ' + '-d%d '*len(days) + '--GWslope %d --tditype %s --tdigen %s --whichtdi %s --lmax %d --f0 %f --df %f --Nf %d %s %s %s ' ) %
#                   tuple( days + [ setup['GWslope'] , setup['tditype'] , setup['tdigen'] , setup['whichtdi'] ,
#                                            setup['lmax'] , setup['f0'] , setup['df'] , setup['Nf'] ,
#                                            setup['Ppath'] , setup['orfdir'] , csddir ] ) )
    print 'done'


