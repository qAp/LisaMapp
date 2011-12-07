#!/usr/bin/env python
import os
import sys
import glob
import re
import cPickle as cpkl
import numpy as np
import myLISAmodule as mlisar
import AnisotropySearch as AS
import optparse


parser = optparse.OptionParser( 'Usage: ./x_simulate_noise.py SETUP.pkl' )
( options , args ) = parser.parse_args()
if len( args ) < 1 :
    parser.error( 'You must specify a SETUP.pkl containing input parameters!' )
file = open( args[0] , 'rb' ) ; setup = cpkl.load( file ) ; file.close()


tditype , tdigen , whichtdis = setup['tditype'] , setup['tdigen'] , setup['whichtdis']
stime , duration , inittime , seed = setup['stime'] , setup['duration'] , setup['inittime'] , setup['seed']
tspath = setup['tspath']


tsdir = '/'.join( re.split( '/' , tspath )[:-1] ) + '/'
print tsdir
if not tspath in glob.glob( tspath ) :
    t , n = mlisar.get_tdiNoise_freq_domain_CovarMatrix( tditype , tdigen , whichtdis ,
                                                         stime , duration , inittime , seed )
    tscale = { 'Cadence1':stime , 'Offset1':inittime }
    tsdict = {}
    tsdict['1'] , tsdict['2'] , tsdict['3'] = AS.Coarsable( n[0] , **tscale ) , AS.Coarsable( n[1] , **tscale ) , AS.Coarsable( n[2] , **tscale )
    tsdict['t'] = AS.Coarsable( t , **tscale )

    if tsdir not in glob.glob( tsdir ) :
        os.system( 'mkdir %s' % tsdir )
    print 'saving time-series to disk'
    file = open( tspath , 'wb' )
    cpkl.dump( tsdict , file , -1 ) ; file.close()
else :
    print '%s already saved to disk.  Nothing to do.' % tspath

