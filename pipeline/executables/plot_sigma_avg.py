#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser

usage = """
%prog AVGSIGDIR AVGSIGFIGDIR\n
AVGSIGDIR --- path to directory containing the average sigma\n
AVGSIGFIGDIR --- path to directory to save the plot of average sigma vs lmax in
"""
parser = OptionParser( usage = usage )

( options , args ) = parser.parse_args()

if len( args ) < 2 :
    parser.error( 'You must specify at least AVGSIGDIR and AVGSIGFIGDIR! For help use -h.' )

avgsigdir , avgsigfigdir = args
avgsigpath = avgsigdir + 'sigma_avg.pkl'
file = open( avgsigpath , 'rb' ) ; sigma_avgs = cpkl.load( file ) ; file.close()

lmaxs = np.arange( len( sigma_avgs ) )
print lmaxs
print sigma_avgs

fig = plt.figure()
ax = fig.add_subplot( 111 )
#ax.plot( lmaxs , sigma_avgs , lmaxs , 1./(lmaxs+1)**2 )
ax.plot( lmaxs , sigma_avgs )
ax.set_title( 'average sigma vs lmax' )
ax.set_ylabel( r'$\sigma_{avg}$' )
ax.set_xlabel( 'lmax' )
if avgsigfigdir not in glob.glob( avgsigfigdir ) :
    os.system( 'mkdir -p %s' % avgsigfigdir )
fig.savefig( avgsigfigdir + 'sigma_avg.png' )

