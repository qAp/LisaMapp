#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import AnisotropySearch as AS

GWslopes = [ 0 ]
IJs = [ 'AE' , 'AT'  , 'ET' , 'AE_AT' , 'AE_ET' , 'AT_ET' , 'AE_AT_ET' ]
norm = 1e-40 
P0 = 1e-34
P00 = P0 * np.sqrt( 4 * np.pi )
print 'expected P0 = %e' % P0
print 'expected P00 = %e' % P00

workdir = os.getcwd() + '/'

for slope in GWslopes :
    for IJ in IJs :
        print 'GWslope = %d , IJ = %s' % ( slope , IJ )
        os.chdir( workdir )
        Ppath = workdir + 'GW_slope_%d/%s/P/P.pkl' % ( slope , IJ )
        stdPpath = workdir + 'GW_slope_%d/%s/stdP/stdP.pkl' % ( slope , IJ )
        sigmamappath = workdir + 'GW_slope_%d/%s/sigmamap/S.pkl' % ( slope , IJ )
        file = open( Ppath , 'rb' ) ; P = cpkl.load( file ) ; file.close()
        file = open( stdPpath , 'rb' ) ; stdP = cpkl.load( file ) ; file.close()
        file = open( sigmamappath , 'rb' ) ; S = cpkl.load( file ) ; file.close()
        print '%20s'*7  % ( 'estimated P00' , 'dP00' , 'SP00' , 'stdP00' , 'estimated P0' , 'dP0' , 'SP0' )
        print '%20e'*7 % ( np.real( P.plm[0] )*norm , ( np.real( P.plm[0] )*norm - P00 ) , np.real( S.plm[0] )*norm , np.real( stdP[0] )*norm ,
                          P.P[ 0,0 ]*norm , ( P.P[ 0,0 ]*norm - P0 ) , S.P[ 0,0 ]*norm )
        print '-'*20*7
