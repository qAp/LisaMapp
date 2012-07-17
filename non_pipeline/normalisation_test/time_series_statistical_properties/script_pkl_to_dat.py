#!/usr/bin/env python
import os
import sys
import cPickle as cpkl

statsname = 'stats.pkl'

file = open( statsname , 'rb' ) ; stats = cpkl.load( file ) ; file.close()

N_params , N_days = stats.shape

file = open( 'stats.dat' , 'w' )
print >> file , '#| DAY | mu_y1 | mu_y2 | sig_y1y1 | sig_y1y2 | sig_y2y1 | sig_y2y2 |'
#print >> file , ( '%15s'*N_params ) % ('mu_y1' ,'mu_y2' ,'sig_y1y1' , 'sig_y1y2' , 'sig_y2y1' , 'sig_y2y2' )
for d in range( N_days ) :
    print >> file , ( '%15d' + '%15e'*N_params ) % tuple( [d+1] + stats[ : , d ].flatten().tolist() )
file.close()

