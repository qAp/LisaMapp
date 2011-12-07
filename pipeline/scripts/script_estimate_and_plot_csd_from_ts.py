#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import myLISAmodule as mlisar
import pylab as pl
import matplotlib.pyplot as plt


days = [ 66 ]

tsdir = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/runs/EccenctricInclined_eta0_0_xi0_0_sw_1_t0_0/signal/point_sources/lon_263_lat_-35_plus_P00_GWslope_0/TDI/Michelson/G2/AET/simulation/stime_1.0_N_86400_for_days/data/'


figdir = 'figures_csd_non_pipeline/'

Nseg = 100


for day in days :
    print 'Working on day %d' % day

    tspath = tsdir + 'd%03d.pkl' % day

    file = open( tspath , 'rb' )
    tsdict = cpkl.load( file ) ; file.close()

    t = tsdict['t'].data
    n1 , n2 , n3 = tsdict['1'].data , tsdict['2'].data , tsdict['3'].data

    stime = t[1] - t[0] ; inittime = t[0] ; N = t.shape[0]

    print "Taking %d averages to estimate spectral densities ... " % ( 2*Nseg )
    nfft = int( np.round( N / Nseg ) )
    noverlap = nfft / 2
    fs = 1. / stime

    P12 , f = pl.csd( n1 , n2 , nfft , fs , noverlap = noverlap )
    P13 , f = pl.csd( n1 , n3 , nfft , fs , noverlap = noverlap )
    P23 , f = pl.csd( n2 , n3 , nfft , fs , noverlap = noverlap )

    fig = plt.figure()
    ax = fig.add_subplot( 311 )
    ax.plot( f , np.real( P12 ) , f , np.imag( P12 ) )
    ax.legend( ( 'real' , 'imag' ) )
    ax.set_ylabel( 'P12' )
    ax.set_title( 'Day %d' % day )
    ax = fig.add_subplot( 312 )
    ax.plot( f , np.real( P13 ) , f , np.imag( P13 ) )
    ax.legend( ( 'real' , 'imag' ) )
    ax.set_ylabel( 'P13' )
    ax = fig.add_subplot( 313 )
    ax.plot( f , np.real( P23 ) , f , np.imag( P23 ) )
    ax.legend( ( 'real' , 'imag' ) )
    ax.set_ylabel( 'P23' ) ; ax.set_xlabel( 'freq[Hz]' )
    if figdir not in glob.glob( figdir ) :
        os.system( 'mkdir -p %s' % figdir )
    fig.savefig( figdir + 'd%03d.png' % day )
    plt.show()
