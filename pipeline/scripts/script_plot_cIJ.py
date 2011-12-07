#!/usr/bin/env python
import os
import sys
import glob 
import cPickle as cpkl
import numpy as np
import matplotlib.pyplot as plt


days = [ 66 ]#, 166 , 266 ]
cIJdir = 'cIJ/'
figdir = 'figures_cIJ/'
norm = 1.





if figdir not in glob.glob( figdir ) :
    os.system( 'mkdir %s' % figdir )

for day in days :
    print "~~~~~~~ Day %03d ~~~~~~~" % day
    file = open( 'cIJ/d%03d.pkl' % day , 'rb' ) ; cIJdict = cpkl.load( file ) ; file.close()

    f = cIJdict['f'].data 
    PII = cIJdict['pII'].data * norm
    PJJ = cIJdict['pJJ'].data * norm
    PIJ = cIJdict['csdata'].data * norm

    fig = plt.figure()
    fig.suptitle( 'cIJs. Day %d' % day )
    ax1 = fig.add_subplot( 311 )
    ax1.plot( f , PII )
    ax1.set_ylabel( 'PII' ) 
    ax2 = fig.add_subplot( 312 )
    ax2.plot( f , PJJ )
    ax2.set_ylabel( 'PJJ' )
    ax3 = fig.add_subplot( 313 )
    ax3.plot( f , np.real( PIJ ) , f , np.imag( PIJ ) )
    ax3.set_ylabel( 'PIJ' ) ; ax3.set_xlabel( 'frequency[Hz]' )
    plt.savefig( figdir + 'd%03d.png' % day )
