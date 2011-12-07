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

figdir = 'figures_psd_non_pipeline/'

Nseg = 50


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

    P11 , f = pl.psd( n1 , nfft , fs , noverlap = noverlap )
    P22 , f = pl.psd( n2 , nfft , fs , noverlap = noverlap )
    P33 , f = pl.psd( n3 , nfft , fs , noverlap = noverlap )

#    print 'Getting theoretical PSD curves...'
#    P11q = mlisar.get_tdiNSD( 'Michelson' , 'G2' , 'A' , 'A' , f )
#    P22q = mlisar.get_tdiNSD( 'Michelson' , 'G2' , 'E' , 'E' , f )
#    P33q = mlisar.get_tdiNSD( 'Michelson' , 'G2' , 'T' , 'T' , f )
#    print "done"

    fig = plt.figure()
    ax = fig.add_subplot( 311 )
    ax.plot( f , P11 )
    ax.set_ylabel( 'P11' )
    ax.set_title( 'Day %d' % day )
    ax = fig.add_subplot( 312 )
    ax.plot( f , P22 )
    ax.set_ylabel( 'P22' )
    ax = fig.add_subplot( 313 )
    ax.plot( f , P33 )
    ax.set_ylabel( 'P33' ) ; ax.set_xlabel( 'freq[Hz]' )
    if figdir not in glob.glob( figdir ) :
        os.system( 'mkdir -p %s' % figdir )
    fig.savefig( figdir + 'd%03d.png' % day )
#    plt.show()
