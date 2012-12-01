#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import matplotlib.pyplot as plt
from optparse import OptionParser
import AnisotropySearch as AS



usage = """
%prog FISHPATH FIGPATH\n
FISHPATH --- path to the file in which Fisher matrix is saved\n
FIGPATH --- path to the file to save the figure
"""
parser = OptionParser( usage = usage )
parser.add_option( '--lmax' , action='store' , dest='lmax' , type='int' , default=15 ,
                   help="l up to which to truncate the original Fisher matrix." )

( options , args ) = parser.parse_args()

if len( args ) < 2 :
    print 'You must specify FISHPATH and FIGPATH!'

fishpath , figpath = args
figdir = os.path.dirname( figpath )

fish = AS.FisherMatrix( fishpath , lmax = options.lmax )
fish.svd()
fish.s

fig = plt.figure()
ax = fig.add_subplot(111)
ax.semilogy( fish.s )
ax.set_ylabel( 'singular value' ) ; ax.set_xlabel( 'number' )
ax.set_title( 'Singular values. lmax = %d' % options.lmax )
if figdir not in glob.glob( figdir ) :
    os.system( 'mkdir -p %s' % figdir )
fig.savefig( figpath )




