#!/usr/bin/env python
import os
import sys
import glob 
import cPickle as cpkl
import matplotlib ; matplotlib.use( 'Agg' )
import matplotlib.pyplot as plt


#sourcename = 'noise'
days = range( 4 , 365 + 1 )
psddir = 'psd/'
figdir = 'figures_psd/'







if figdir not in glob.glob( figdir ) :
    os.system( 'mkdir %s' % figdir )

for day in days :

    print "~~~~~~~ Day %03d ~~~~~~~" % day

    psdpath = psddir + 'd%03d.pkl' % day
    file = open( psdpath , 'rb' ) ; psddict = cpkl.load( file ) ; file.close()

    f = psddict['f'].data
    P11 = psddict['AA'].data
    P22 = psddict['EE'].data
    P33 = psddict['TT'].data

    figname = 'd%03d.png' % day

    fig = plt.figure()
    ax1 = fig.add_subplot( 311 )
    ax1.plot( f , P11 )
    ax1.set_ylabel( 'P11' )
    ax1.set_title( 'Power spectral density' )
    ax2 = fig.add_subplot( 312 )
    ax2.plot( f , P22 )
    ax2.set_ylabel( 'P22' )
    ax3 = fig.add_subplot( 313 )
    ax3.plot( f , P33 )
    ax3.set_ylabel( 'P33' )
    ax3.set_xlabel( 'frequency $[Hz]$' )
    fig.savefig( figdir + figname )
