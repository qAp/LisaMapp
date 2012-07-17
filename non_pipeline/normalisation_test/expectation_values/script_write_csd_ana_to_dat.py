#!/usr/bin/env python
import os
import sys
import cPickle as cpkl
import numpy as np

file = open( 'csd_ana/d001.pkl' , 'rb' )
csddict = cpkl.load( file ) ; file.close()

Nf = csddict['f'].data.shape[0]

file = open( 'csd_ana.dat' , 'w' )
print >> file , ( '#' + '%19s' + 2*'%20s' ) % ( 'frequency' , 'analytical Re{csd}' , 'analytical Im{csd}' )
for k in range( Nf ) :
    print >> file , ( 3*'%20f' ) % ( csddict['f'].data[k] ,
                                     np.real(csddict['AE'].data[k]) , np.imag(csddict['AE'].data[k]) )
file.close()
