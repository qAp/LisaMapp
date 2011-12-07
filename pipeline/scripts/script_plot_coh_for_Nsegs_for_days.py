
#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import myLISAmodule as mlisar
import pylab as pl
import matplotlib.pyplot as plt


days = range( 1 , 365+1 ) ; days.remove( 1 ) ; days.remove( 180 ) ; days.remove( 270 )

tsdir = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/runs/EccenctricInclined_eta0_0_xi0_0_sw_1_t0_0/noise/TDI/Michelson/G2/AET/simulation/by_me/stime_1.0_N_86400_for_days/data/'


for day in days :
    print 'Working on day %d' % day

    tspath = tsdir + 'd%03d.pkl' % day

    file = open( tspath , 'rb' )
    tsdict = cpkl.load( file ) ; file.close()

    t = tsdict['t'].data
    n1 , n2 , n3 = tsdict['s1'].data , tsdict['s2'].data , tsdict['s3'].data

    stime = t[1] - t[0] ; inittime = t[0] ; N = t.shape[0]
    fs = 1. / stime

    power_max = 3
    Nsegs = 10**np.arange( power_max + 1 )
    P12s = [] ; P12s_legend = [] ; C12s = [] ; C12s_legend = []
    P13s = [] ; P13s_legend = [] ; C13s = [] ; C13s_legend = []
    P23s = [] ; P23s_legend = [] ; C23s = [] ; C23s_legend = []
    for i , Nseg in enumerate( Nsegs ) :
        print "Taking %d averages to esimate spectral densities ... " % ( 2*Nseg )
        nfft = int( np.round( N / Nseg ) )
        noverlap = nfft / 2

        P11 , f = pl.psd( n1 , nfft , fs , noverlap = noverlap )
        P22 , f = pl.psd( n2 , nfft , fs , noverlap = noverlap )
        P33 , f = pl.psd( n3 , nfft , fs , noverlap = noverlap )

        P12 , f = pl.csd( n1 , n2 , nfft , fs , noverlap = noverlap )
        P23 , f = pl.csd( n2 , n3 , nfft , fs , noverlap = noverlap )
        P13 , f = pl.csd( n1 , n3 , nfft , fs , noverlap = noverlap )

        C12 = np.abs( P12 )**2 / ( P11*P22 )
        C23 = np.abs( P23 )**2 / ( P22*P33 )
        C13 = np.abs( P13 )**2 / ( P11*P33 )
        for k in range( C12.shape[0] ) :
            if C12[k] < 0 :
                print 'C12[%d] < 0' % k
        P12s += [ f , np.abs( P12 ) ] ; C12s += [ f , C12 ] ; C12s_legend += [ '%d' % Nseg ]
        P23s += [ f , np.abs( P23 ) ] ; C23s += [ f , C23 ] ; C23s_legend += [ '%d' % Nseg ]
        P13s += [ f , np.abs( P13 ) ] ; C13s += [ f , C13 ] ; C13s_legend += [ '%d' % Nseg ]
        print 'done'
        
    print 'C12'
    print C12s[0].shape , C12s[1].shape
    print C12s[2].shape , C12s[3].shape
    print C12s[4].shape , C12s[5].shape

    figdir = 'figures_coherences/day%03d/' % day
    if figdir not in glob.glob( figdir ) :
        os.system( 'mkdir -p %s' % figdir )

    fig = plt.figure()
    ax = fig.add_subplot( 111 )
#    ax.semilogy( *tuple( C12s ) )
    h1 = ax.semilogy( C12s[0] , C12s[1] )
    ax.hold( True )
    h2 = ax.semilogy( C12s[2] , C12s[3] )
    ax.hold( True )
    h3 = ax.semilogy( C12s[4] , C12s[5] )
    ax.hold( True )
    h4 = ax.semilogy( C12s[6] , C12s[7] )
    ax.hold( False )
    ax.legend( ( h1 , h2 , h3 , h4 ) , C12s_legend )
    ax.set_ylabel( 'C12' ) ; ax.set_xlabel( 'freq[Hz]' )
    ax.set_title( 'Day %d' % day )
    fig.savefig( figdir + 'coh12_d%03d.png' % day )

    fig = plt.figure()
    ax = fig.add_subplot( 111 )
#    ax.semilogy( *tuple( C13s ) )
    h1 = ax.semilogy( C13s[0] , C13s[1] )
    ax.hold( True )
    h2 = ax.semilogy( C13s[2] , C13s[3] )
    ax.hold( True )
    h3 = ax.semilogy( C13s[4] , C13s[5] )
    ax.hold( True )
    h4 = ax.semilogy( C13s[6] , C13s[7] )
    ax.hold( False )
    ax.legend( ( h1 , h2 , h3 , h4 ) , C13s_legend )
    ax.set_ylabel( 'C13' ) ; ax.set_xlabel( 'freq[Hz]' )
    ax.set_title( 'Day %d' % day )
    fig.savefig( figdir + 'coh13_d%03d.png' % day )    

    fig = plt.figure()
    ax = fig.add_subplot( 111 )
#    ax.semilogy( *tuple( C23s ) )
    h1 = ax.semilogy( C23s[0] , C23s[1] )
    ax.hold( True )
    h2 = ax.semilogy( C23s[2] , C23s[3] )
    ax.hold( True )
    h3 = ax.semilogy( C23s[4] , C23s[5] )
    ax.hold( True )
    h4 = ax.semilogy( C23s[6] , C23s[7] )
    ax.hold( False )
    ax.legend( ( h1 , h2 , h3 , h4 ) , C23s_legend )
    ax.set_ylabel( 'C23' ) ; ax.set_xlabel( 'freq[Hz]' )
    ax.set_title( 'Day %d' % day )
    fig.savefig( figdir + 'coh23_d%03d.png' % day )





#        
#    fig12 = plt.figure()
#    ax12 = fig12.add_subplot( 111 )
#    ax12.semilogy( *tuple( C12s ) ) ; ax12.legend( tuple( C12s_legend ) )
#    ax12.set_ylabel( 'C12' ) ; ax12.set_xlabel( 'frequency[Hz]' ) #; ax12.grid( color='c' )
#    fig12.savefig( figdir + 'coh12_d%03d.png' % day )
#    
#    fig13 = plt.figure()
#    ax = fig13.add_subplot( 111 )
#    ax.semilogy( *C13s ) ; ax.legend( C13s_legend )
#    ax.set_ylabel( 'C13' ) ; ax.set_xlabel( 'frequency[Hz]' ) ; ax.grid( color='c' )
#    fig13.savefig( figdir + 'coh13_d%03d.png' % day )
#    
#    fig23 = plt.figure()
#    ax = fig23.add_subplot( 111 )
#    ax.semilogy( *C23s ) ; ax.legend( C23s_legend )
#    ax.set_ylabel( 'C23' ) ; ax.set_xlabel( 'frequency[Hz]' ) ; ax.grid( color='c' )
#    fig23.savefig( figdir + 'coh23_%03d.png' % day )






    
