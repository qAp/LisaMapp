#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
import AnisotropySearch as AS
import PostProcess as PP


usage = """
%prog STDPPATH FIGPATH\n
STDPPATH --- path to file containing standard deviations of Plm\n
FIGPATH --- path to file to save the figure in
"""
parser = OptionParser()

parser.add_option( '--norm' , action='store' , dest='norm' , type='string' , default='1.' ,
                   help='Factor by which the standard deviations should bey multiplied before being plotted.' )

( options , args ) = parser.parse_args()

if len( args ) < 2 :
    parser.error( 'You must specify STDPPATH and FIGPATH! ' )

stdPpath , figpath = args[ :2 ]

file = open( stdPpath , 'rb' ) ; stdP = cpkl.load( file ) ; file.close()

lmax = int( np.sqrt( len( stdP ) ) ) - 1

dmls = [ 2*l+1 for l in range( lmax+1 ) ]
xticksat = np.array( [ sum(dmls[:ii]) for ii in range(len(dmls)) ] )
xticksare = np.array( range( lmax+1 ) )

plt.figure()
ax = plt.subplot(111)
ax.semilogy( float( options.norm ) * stdP )
plt.xticks( xticksat , xticksare , rotation=45 )
ax.grid( b='on' )
plt.xlim( ( 0 , stdP.shape[0] ) )
plt.ylabel( 'standard deviation of $P_{lm}$' )
plt.xlabel( 'l ((l,m) packed as $ l^{2} + l + m + 1 $)' )
figdir = os.path.dirname( figpath )
if figdir not in glob.glob( figdir ) :
    os.system( 'mkdir -p %s' % figdir )
plt.savefig( figpath )
    
