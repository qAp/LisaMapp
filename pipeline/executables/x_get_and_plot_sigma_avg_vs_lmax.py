#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
import AnisotropySearch as AS



parser = OptionParser( "usage: ./x_get_and_plot_sigma_avg_vs_lmax.py SETUP.pkl" )
( options , args ) = parser.parse_args()
if len( args ) < 1 :
    parser.error( 'You must specify a SETUP file of parameters!' )
setupname = args[0]
file = open( setupname , 'rb' ) ; setup = cpkl.load( file ) ; file.close()

lmax_max = setup['maxlike']['lmax']
Gdir = setup['maxlike']['Gdir']
figdir = setup['sigma_avg']['figdir']
nlon , nlat = setup['sigma_avg']['nlon'] , setup['sigma_avg']['nlat']



lmaxs = range( lmax_max + 1 )
sigma_avgs = []
for lmax in lmaxs :
    print 'Calculating average sigma at lmax = %d' % lmax
    
    fish = AS.FisherMatrix( Gdir + 'G.pkl' , lmax = lmax )
    fishinv = fish.invert()
    map_sigma , lats , lons = AS.getSigmaMap( fishinv , nlat , nlon )

    dlon = lons[1] - lons[0] ; dlat = lats[0] - lats[nlon]
    dOmega = np.radians( dlon ) * np.radians( dlat ) * np.sin( np.pi/2 - np.radians( lats ) )

    sigma_avg = np.sum( dOmega * map_sigma ) / ( lmax + 1 )**2
    sigma_avgs += [ sigma_avg ]


if figdir not in glob.glob( figdir ) :
    os.system( 'mkdir %s' % figdir )
    
fig = plt.figure()
ax = fig.add_subplot( 111 )
ax.plot( lmaxs , sigma_avgs )
ax.set_title( 'average standard deviation vs lmax' )
ax.set_ylabel( r'$\sigma_{avg}$' ) ; ax.set_xlabel( 'lmax' )
fig.savefig( figdir + 'sigma_avg_vs_lmax.png' )

