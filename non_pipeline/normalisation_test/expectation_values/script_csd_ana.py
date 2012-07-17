#!/usr/bin/env python
import os
import sys
import glob
import numpy as np
import cPickle as cpkl
import synthlisa
import myLISAmodule as mlisar
import AnisotropySearch as AS


days = range( 1 , 365+1 )
csddir = 'csd_ana'
P00 = 1/90.
stime = 8.
df = 1 / 5760.
tditype , tdigen = 'Michelson' , 'G2'




fs = 1 / stime
f0 , Nf = 0 , int( (fs/2) / df ) + 1
f = f0 + df * np.arange( Nf )
fscale = { 'Offset1':f0 , 'Cadence1':df }


csddict = {}
IJs = [ 'AE' , 'AT' , 'ET' ]
for IJ in IJs :
    g00 = np.ones( f.shape ) * 288.
    pIJ = g00 * P00
    csddict[IJ] = AS.Coarsable( np.copy( pIJ ) , **fscale )

csddict['f'] = AS.Coarsable( np.copy( f ) , **fscale )


if csddir not in glob.glob( csddir ) :
    os.system( 'mkdir %s' % csddir )
    
for day in days :
    file = open( csddir + '/d%03d.pkl' % day , 'wb' ) ; cpkl.dump( csddict , file , -1 ) ; file.close()


