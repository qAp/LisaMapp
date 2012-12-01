#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import optparse
import matplotlib.pyplot as plt


usage = """
%prog TSPATH FIGPATH\n
TSPATH --- path to file for the time-series\n
FIGPATH --- path to file to save the figure in
"""
parser = optparse.OptionParser( usage = usage )

parser.add_option( '--scale_ts' , action='store' , dest='scale_ts' , type='float' , default=1. ,
                   help='Scale the time-series by this number before plotting' )

( options , args ) = parser.parse_args()

if len( args ) < 2 :
    parser.error( "You must specify at least TSPATH and FIGPATH!  Type './plot_ts.py -h'" )
else :
    tspath , figpath = args[ :2 ]

file = open( tspath , 'rb' ) ; tsdict = cpkl.load( file ) ; file.close()

t = tsdict['t'].data
s1 = options.scale_ts * tsdict['1'].data
s2 = options.scale_ts * tsdict['2'].data
s3 = options.scale_ts * tsdict['3'].data

fig = plt.figure()
fig.suptitle( 'Time-Series' )
ax1 = fig.add_subplot( 311 )
ax1.plot( t , s1 )
ax1.set_ylabel( 's1' )
ax2 = fig.add_subplot( 312 )
ax2.plot( t , s2 )
ax2.set_ylabel( 's2' )
ax3 = fig.add_subplot( 313 )
ax3.plot( t , s3 )
ax3.set_ylabel( 's3' ) ; ax3.set_xlabel( 'time[s]' )
figdir = os.path.dirname( figpath )
if figdir not in glob.glob( figdir ) :
    os.system( 'mkdir -p %s' % figdir )
fig.savefig( figpath )


