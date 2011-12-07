#!/usr/bin/env python
import os
import sys
import glob 
import cPickle as cpkl
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt


days = [ 66 , 166 , 266 ]
psddir = 'psd/'
figdir = 'figures_psd/'



if figdir not in glob.glob( figdir ) :
    os.system( 'mkdir %s' % figdir )

for day in days :
    print "~~~~~~~ Day %03d ~~~~~~~" % day
    
    file = open( psddir + 'd%03d.pkl' % day , 'rb' ) ; psddict = cpkl.load( file ) ; file.close()
    f = psddict['f'].data
    P11 = psddict['AA'].data
    P22 = psddict['EE'].data
    P33 = psddict['TT'].data

    print 'Estimating power spectral densities from some simulated time-series...'
    tsdirz = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/runs/EccenctricInclined_eta0_0_xi0_0_sw_1_t0_0/output/signal/point_sources/lon_263_lat_-35_plus_P00_GWslope_0/TDI/Michelson/G2/AET/noise/TDI/Michelson/G2/AET/simulation/by_me/stime_1.0_N_86400_for_days/data/'
    tspathz = tsdirz + 'd%03d.pkl' % day
    file = open( tspathz , 'rb' ) ; tsdictz = cpkl.load( file ) ; file.close()
    tz = tsdictz['t'].data
    s1z , s2z , s3z = tsdictz['1'].data , tsdictz['2'].data , tsdictz['3'].data
    stimez = tz[1] - tz[0] ; inittimez = tz[0] ; Nz = tz.shape[0]
    nfftz = int( np.round( Nz / 100 ) )
    noverlapz = nfftz / 2
    fsz = 1. / stimez
    P11z , fz = pl.psd( s1z , nfftz , fsz , noverlap = noverlapz )
    P22z , fz = pl.psd( s2z , nfftz , fsz , noverlap = noverlapz )
    P33z , fz = pl.psd( s3z , nfftz , fsz , noverlap = noverlapz )
    print 'done'

    fig = plt.figure()
    ax1 = fig.add_subplot( 311 )
    h1 = ax1.plot( fz , P11z , f , P11 )
#    h1 = ax1.plot( f , P11 )
    ax1.set_ylabel( 'P11' )
    ax1.set_title( 'Power Spectral Densities' )
    ax2 = fig.add_subplot( 312 )
    h2 = ax2.plot( fz , P22z , f , P22 )
#    h2 = ax2.plot( f , P22 )
    ax2.set_ylabel( 'P22' )
    ax3 = fig.add_subplot( 313 )
    h3 = plt.plot( fz , P33z , f , P33 )
#    h3 = plt.plot( f , P33 )    
    ax3.set_ylabel( 'P33' )
    ax3.set_xlabel( 'frequency$[Hz]$' )
    plt.savefig( figdir + 'd%03d.png' % day )
