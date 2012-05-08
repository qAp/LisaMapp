#!/usr/bin/env python
import os
import sys
import glob
import optparse
import cPickle as cpkl
import numpy as np
import AnisotropySearch as AS
import myUsefuls as mufls


def get_comatrix( f ) :
    """ Define the covariance matrix here!!! """
    Nf = f.shape[ 0 ] 
    p12 , p13 , p23 = np.zeros( (Nf,) ) , np.zeros( (Nf,) ) , np.zeros( (Nf,) )
    p11 , p22 , p33 = 120. * np.ones( (Nf,) ) , 120. * np.ones( (Nf,) ) , 120. * np.ones( (Nf,) )
    comatrix = np.array( [ [ p11 , p12 , p13 ] ,
                           [ np.conj(p12) , p22 , p23 ] ,
                           [ np.conj(p13) , np.conj(p23) , p33 ] ] )
    return comatrix


usage = """
%prog TSDIR\n
TSDIR --- directory in which to save the time-series files
"""
parser = optparse.OptionParser( usage=usage )
parser.add_option( '-d' , '--day' , action='append' , dest='days' , type='int' , nargs=1 , 
                   help='Days for which to simulate noise time-series' )
parser.add_option( '--Tseg' , action='store' , dest='Tseg' , type='float' , nargs=1 ,
                   help='Duration of stitched segments making up the output time-series' )
parser.add_option( '--stime' , action='store' , dest='stime' , type='float' , nargs=1 , default='1.0' ,
                   help='sampling time of time-series' )
parser.add_option( '--seed' , action='store' , dest='seed' , type='string' , nargs=1 , default='random' ,
                   help='seed for random number generator' )

( options , args ) = parser.parse_args()

if len( args ) < 1 :
    parser.error( 'You must specify TSDIR!' )
else :
    tsdir = args[ 0 ]

dayinsecs = 86400. 
Nvar = 3

N = int( dayinsecs / options.stime )
n = int( options.Tseg / options.stime )

if N % n == 0 :
    Nseg = int( N / n ) * 2
else :
    Nseg = ( int( N / n ) + 1 ) * 2

duration = n * options.stime

freqdict = mufls.get_freqs_from_duration_and_stime( options.stime , duration )

comatrix = get_comatrix( freqdict['f'] )

for day in options.days :
    print 'Day %d' % day
    t0 = ( day - 1 )*dayinsecs
    t = t0 + options.stime * np.arange( N )
    Npd_before_today = (day - 1) * Nseg * 2 * Nvar * freqdict['Nf'] 

    s = 1 ; tails = [ None ]*Nvar ; TSs = [ [] for v in range( Nvar ) ]
    while s <= Nseg :
        print 'Segments (%d|%d)' % ( s , s+1 )
        Npd_before_segl = Npd_before_today + (s-1)*2*Nvar*freqdict['Nf']
        Npd_before_segr = Npd_before_today + s*2*Nvar*freqdict['Nf']
        tl , tsl = mufls.get_noise_freq_domain_CovarMatrix(comatrix,freqdict['df'],
                                                           t0,freqdict['parityN'],
                                                           options.seed,
                                                           Npd_before_segl )
        tr , tsr = mufls.get_noise_freq_domain_CovarMatrix(comatrix,freqdict['df'],
                                                           t0,freqdict['parityN'],
                                                           options.seed,
                                                           Npd_before_segr )

        for v in range( Nvar ) :
            print 'v = ' , v
            if s == 1 :
                print 'tails[v]' , tails[v]
            else :
                print 'tails[v][:3]' , tails[v][:3]
            ts , tail = mufls.window_and_join( tsl[v] , tsr[v] , tails[v] )
            TSs[ v ] +=  list( ts )  ; tails[ v ] = np.copy( tail )
            print 'tails[v][:3]' , tails[v][:3]
        s += 2

    TSs = np.array( TSs )[ : , :N ]
    
    tscale = { 'Cadence1':options.stime , 'Offset1':t0 }
    tsdict = { 't':AS.Coarsable( t , **tscale ) , 
               '1':AS.Coarsable( TSs[0] , **tscale ) ,
               '2':AS.Coarsable( TSs[1] , **tscale ) ,
               '3':AS.Coarsable( TSs[2] , **tscale ) }    
    if tsdir == '' :
        pass
    elif tsdir not in glob.glob( tsdir ) :
        os.system( 'mkdir -p %s' % tsdir )
    print 'saving time-series to disk...'
    file = open( tsdir+'/d%03d.pkl' % day , 'wb' )
    cpkl.dump( tsdict , file , -1 ) ; file.close()    


