#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import AnisotropySearch as AS
import myLISAmodule as mlisar

""" """
stime = .5
df = 1. / 5760. #2880 
lmax = 0
g00 = 9.
IJ = 'EE'
days = range( 1 , 365+1 )
fs = 1 / stime ; fNyq = fs / 2
f0 = df 
Nf = int( ( fNyq - f0 ) / df ) 
""" """



f = f0 + df * np.arange( Nf )

indxp = AS.getMLvec( lmax , 'p' )
SpHreal = np.ones( ( len( indxp ) , Nf ) ) * ( g00/np.sqrt(2*np.pi) )
SpHimag = np.ones( ( len( indxp ) , Nf ) ) * 0.

orfdict = {'OrfMultipleMoments': { 'ntrunc': lmax , 'f': f , 'Antenna': IJ , 'real': SpHreal , 'imag': SpHimag } }


for day in days :
    orfpath = 'tdiI_%s_tdiJ_%s_lmax_%d_f0_%.6f_df_%.6f_Nf_%d_g00_9/data/orf_d%03d.pkl' % ( IJ[0] , IJ[1] , lmax , f0 , df , Nf , day )
    orfdir = os.path.dirname( orfpath )
    if orfdir not in glob.glob( orfdir ) :
        os.system( 'mkdir -p %s' % orfdir )
    file = open( orfpath , 'wb' ) ; cpkl.dump( orfdict , file , -1 ) ; file.close()
    print " Orf multipole moments saved in %s " % orfpath



