#!/usr/bin/env python
import os
import sys
import glob 
import cPickle as cpkl
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt


days = [ 66 , 166 , 266 ]
csddir = 'csd/'
figdir = 'figures_csd/'



if figdir not in glob.glob( figdir ) :
    os.system( 'mkdir %s' % figdir )

for day in days :
    print "~~~~~~~ Day %03d ~~~~~~~" % day

    file = open( csddir + 'd%03d.pkl' % day , 'rb' ) ; csddict = cpkl.load( file ) ; file.close()
    f = csddict['f'].data
    P12 = csddict['AE'].data
    P13 = csddict['AT'].data
    P23 = csddict['ET'].data

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
    P12z , fz = pl.csd( s1z , s2z , nfftz , fsz , noverlap=noverlapz )
    P13z , fz = pl.csd( s1z , s3z , nfftz , fsz , noverlap=noverlapz )
    P23z , fz = pl.csd( s2z , s3z , nfftz , fsz , noverlap=noverlapz )
    print 'done'

    fig = plt.figure()
    fig.suptitle( 'Cross Spectral Densities' )
    ax1 = fig.add_subplot(311)
    ax1.plot( fz , np.real( P12z ) , fz , np.imag( P12z ) , f , np.real( P12 ) , f , np.imag( P12 ) , 'y' )
#    ax1.plot( f , np.real( P12 ) , f , np.imag( P12 ) , 'y' )
#    ax1.plot( fz , np.real( P12z ) , fz , np.imag( P12z ) )
    ax1.legend( ( 'real' , 'imag' ) )
    ax1.set_ylabel( 'P12' )
    ax2 = fig.add_subplot(312)
    ax2.plot( fz , np.real( P13z ) , fz , np.imag( P13z ) , f , np.real( P13 ) , f , np.imag( P13 ) , 'y' )
#    ax2.plot( f , np.real( P13 ) , f , np.imag( P13 ) , 'y' )
#    ax2.plot( fz , np.real( P13z ) , fz , np.imag( P13z ) )
    ax2.legend( ( 'real' , 'imag' ) )
    ax2.set_ylabel( 'P13' )
    ax3 = fig.add_subplot(313)
    ax3.plot( fz , np.real( P23z ) , fz , np.imag( P23z ) , f , np.real( P23 ) , f , np.imag( P23 ) , 'y' )
#    ax3.plot( f , np.real( P23 ) , f , np.imag( P23 ) , 'y' )
#    ax3.plot( fz , np.real( P23z ) , fz , np.imag( P23z ) )
    ax3.legend( ( 'real' , 'imag' ) )
    ax3.set_ylabel( 'P23' )
    ax3.set_xlabel( 'frequency [Hz]' )
    plt.savefig( figdir + 'd%03d.png' % day )
