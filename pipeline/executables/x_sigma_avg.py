#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
from optparse import OptionParser
import numpy as np
import AnisotropySearch as AS


usage="""
$prog GDIR AVGSIGDIR\n
GDIR --- directory containing Fisher matrix\n
AVGSIGDIR --- directory to save average sigma\n
"""
parser = OptionParser( usage = usage )

parser.add_option( '--regMethod' , action='store' , dest='regMethod' , type='int' , default=1 ,
                   help='The method to regularise the Fisher matrix with.' )

parser.add_option( '--regCutoff' , action='store' , dest='regCutoff' , type='string' , default='9.5e-4' ,
                   help='The cut-off for the regularisation of the Fisher matrix.' )

parser.add_option( '--nlon' , action='store' , dest='nlon' , type='int' , default=180 ,
                    help='Number of longitudes on the sky over which to average the sigma.' )

parser.add_option( '--nlat' , action='store' , dest='nlat' , type='int' , default=91 ,
                    help='Number of latitudes on the sky over which to average the sigma.' )

parser.add_option( '--lmax' , action='store' , dest='lmax' , type='int' , default=20 ,
                    help='l up to which average sigma will be calculated for.' )

( options , args ) = parser.parse_args()

if len( args ) < 2 :
    parser.error( 'You must specify GDIR and AVGSIGDIR! See Help: ./x_sigma_avg.py -h' )

Gdir , avgsigdir = args[ :2 ]

avgsigpath = avgsigdir + 'sigma_avg.pkl'

sigma_avgs = []
for lmax in range( options.lmax + 1 ) :
    print 'Calculating average sigma from Fisher matrix truncated at lmax = %d' % lmax

    fish = AS.FisherMatrix( Gdir + 'G.pkl' , lmax=lmax )
    if options.regMethod == 0 :
        print 'Calculating unregularised inverse of Fisher matrix...'
        fishinv = fish.invert()
        print 'done'
    else :
        print "Calculating regularised inverse of Fisher matrix..."
        fish.regularise( regMethod=options.regMethod , regCutoff=float( options.regCutoff ) )
        fishinv = fish.reginvert()
    
    covarm = np.dot( fishinv , np.dot( fish.fish , fishinv ) )
    map_sigma , lats , lons = AS.getSigmaMap( covarm , options.nlat , options.nlon )

    dlon = lons[ 1 ] - lons[ 0 ] ; dlat = lats[ 0 ] - lats[ options.nlon ]
    dOmega = np.radians( dlon ) * np.radians( dlat ) * np.sin( np.pi/2 - np.radians( lats ) )
    sigma_avg = np.sum( dOmega * map_sigma ) / ( lmax + 1 )**2

    sigma_avgs += [ sigma_avg ]

if avgsigdir not in glob.glob( avgsigdir ) :
    os.system( 'mkdir -p %s' % avgsigdir )
file = open( avgsigpath , 'wb' ) ; cpkl.dump( sigma_avgs , file , -1 ) ; file.close()

    


