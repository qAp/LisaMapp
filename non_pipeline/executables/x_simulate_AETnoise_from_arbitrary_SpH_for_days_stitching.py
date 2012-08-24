#!/usr/bin/env python
import os
import sys
import glob
import optparse
import cPickle as cpkl
import numpy as np
import AnisotropySearch as AS
import myUsefuls as mufls
import testsim as tsim


usage = """
%prog ORFDIR PPATH TSDIR\n
ORFDIR --- directory in which to save, or to look for tdiORF_SpHs\n
PPATH --- file path to the SkyMap object describing the P_{lm}s\n
TSDIR --- directory in which to save the time-series files
"""
parser = optparse.OptionParser( usage=usage )
parser.add_option( '-d' , '--day' , action='append' , dest='days' , type='int' , nargs=1 , 
                   help='Days for which to simulate noise time-series' )
parser.add_option( '--Tseg' , action='store' , dest='Tseg' , type='float' , nargs=1 ,
                   help='Duration of stitched segments making up the output time-series' )
parser.add_option( '--stime' , action='store' , dest='stime' , type='float' , nargs=1 , default='1.0' ,
                   help='sampling time of time-series' )
parser.add_option( '--GWSpectralSlope' , action='store' , dest='GWSpectralSlope' , type='int' , nargs=1 , default=0 ,
                   help='GW background spectral slope' )
parser.add_option( '--H0' , action='store' , dest='H0' , type='float' , nargs=1 , default=1.0 ,
                   help='GW background spectrum normalisation constant' )
parser.add_option( '--lmax' , action='store' , dest='lmax' , type='int' , nargs=1 , default=0 ,
                   help='Maximum degree l in SpH considered' )
parser.add_option( '--seed' , action='store' , dest='seed' , type='string' , nargs=1 , default='random' ,
                   help='seed for random number generator' )
parser.add_option( '--compute_ORF_SpHs' , action='store' , dest='compute_ORF_SpHs' , type='int' , nargs=1 , default=0 ,
                   help='1 (to compute tdiORF_SpHs) or 0 (to load from disk)' )

( options , args ) = parser.parse_args()

if len( args ) < 3 :
    parser.error( 'You must specify ORFDIR, PPATH and TSDIR!' )
else :
    orfdir , Ppath , tsdir = args[ :3 ]

IJs = [ 'AA' , 'AE' , 'AT' , 'EE' , 'ET' , 'TT' ]
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

for day in options.days :
    print 'Day %d' % day
    t0 = ( day - 1 )*dayinsecs
    t = t0 + options.stime * np.arange( N )
    orfpaths = [ orfdir + '/%s/d%03d.pkl' % ( IJ , day ) for IJ in IJs ]
    Npd_before_today = (day - 1) * Nseg * 2 * Nvar * freqdict['Nf'] 

    s = 1 ; tails = [ None ]*Nvar ; TSs = [ [] for v in range( Nvar ) ]
    while s <= Nseg :
        print 'Segments (%d|%d)' % ( s , s+1 )
        Npd_before_segl = Npd_before_today + (s-1)*2*Nvar*freqdict['Nf']
        Npd_before_segr = Npd_before_today + s*2*Nvar*freqdict['Nf']
        tl , tsl = tsim.simulate_AETnoise_from_arbitrary_SpH(duration,options.stime,
                                                             t0 + (s-1)*duration/2,
                                                             options.GWSpectralSlope, options.H0,
                                                             Ppath, options.lmax,
                                                             options.seed,
                                                             Npd_before_segl,
                                                             Nvar,
                                                             options.compute_ORF_SpHs,
                                                             *orfpaths )
        tr , tsr = tsim.simulate_AETnoise_from_arbitrary_SpH(duration,options.stime,
                                                             t0 + s*duration/2,
                                                             options.GWSpectralSlope, options.H0,
                                                             Ppath, options.lmax,
                                                             options.seed,
                                                             Npd_before_segr,
                                                             Nvar,
                                                             options.compute_ORF_SpHs,
                                                             *orfpaths )
        for v in range( Nvar ) :
            ts , tail = mufls.window_and_join( tsl[v] , tsr[v] , tails[v] )
            TSs[ v ] +=  list( ts )  ; tails[ v ] = np.copy( tail )
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
        
