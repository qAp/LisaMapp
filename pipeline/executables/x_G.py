#!/usr/bin/env python
import os
import sys
import glob
import AnisotropySearch as AS
import numpy as np
import cPickle as cpkl
from optparse import OptionParser

usage = """
%prog ORFDIR PSDDIR GDIR\n
TSDIR --- path to directory containing time-series\n
ORFDIR --- path to directory containing ORF multipole moments\n
PSDDIR --- path to directory containing PSDs\n
Gdir --- path to directory to save G in
cIIdir --- path to directory to save cIIs in
"""
parser = OptionParser( usage = usage )

parser.add_option( '--scale_ts' , action='store' , dest='scale_ts' , type='float' , default=1 , nargs=1 ,
                   help='Scale the time-series beforehand.' )

parser.add_option( '-d' , '--day' , action='append' , dest='days' , type='int' , nargs=1 ,
                   help='Days to sum over for maximum likelihood estimator of G' )

parser.add_option( '--GWslope' , action='store' , dest='GWslope' , type='int' , default=0 , nargs=1 ,
                   help='Gravitational spectral slope' )

parser.add_option( '--flow' , action='store' , dest='flow' , type='float' , default=1e-4 , nargs=1 ,
                   help='Lower limit of freqeuncy range for maximum likelihood estimator of G' )

parser.add_option( '--fhigh' , action='store' , dest='fhigh' , type='float' , default=1e-1 , nargs=1 ,
                   help='Higher limit of frequency range for maximum likelihood estimator of G' )

parser.add_option( '--lmax' , action='store' , dest='lmax' , type='int' , default=20 , nargs=1 ,
                   help='Maximum l to take into account in estimating G' )

( options , args ) = parser.parse_args()


if len( args ) < 5 :
    parser.error( 'You must specify TSDIR, ORFDIR, PSDDIR, GDIR and cIIdir! See help ./G.py -h' )

tsdir , orfdir , psddir , Gdir , cIIdir = args[ :5 ] 

Gpath = Gdir + '/G.pkl'

if cIIdir not in glob.glob( cIIdir ) :
    os.system( 'mkdir -p %s' % cIIdir )

firstavailable = True ; days_skipped = []
for day in options.days :
    print "\n~~~ DAY %d ~~~" % day

    tspath  = tsdir  + '/d%03d.pkl' % day
    psdpath = psddir + '/d%03d.pkl' % day
    orfpath = orfdir + '/orf_d%03d.pkl' % day

    if tspath not in glob.glob( tspath ) :
        print 'Time-series not available on day %03d' % day ; days_skipped += [ day ] ; continue
    if psdpath not in glob.glob( psdpath ) :
        print 'PSD not available on day %03d' % day ; days_skipped += [ day ] ; continue
    if orfpath not in glob.glob( orfpath ) :
        print 'ORF not availalbe on day %03d' % day ; days_skipped += [ day ] ; continue

    print "Time-Series, psd and orf all available on day %03d." % day
    print "First available day?" , firstavailable

    file = open( psdpath , 'rb' ) ; psddict = cpkl.load( file ) ; file.close()

    orf = AS.OrfMultipleMoments( orfpath )

    ggg = AS.GGG( orf , psddict , options.GWslope , day , cIIdir=cIIdir  )
    GG = ggg.getsummand( flow = options.flow , fhigh = options.fhigh , lmax= options.lmax )

    if firstavailable:
        G = np.zeros( GG.data.shape , complex )
        print 'Calculating normalisation factor due to coarsegraining and windowing with a Hanning window for both time-series'
        file = open( tspath ) ; tsdict = cpkl.load( file ) ; file.close()
        ts = AS.TimeSeries( tsdict )
        N = ts.t.data.shape[0] ; T = ts.t.Cadence1 * N ; df = ggg.fcoarse.Cadence1
        window = np.hanning( N )
        norm = ( np.sum( window**2 ) / N ) / ( np.sum( window**4 ) / N ) * T*df
        firstavailable = False
        
    G += GG.data

print 'normalise for coarsegraining and windowing with a Hanning window for both time-series'
G = AS.Coarsable( norm * G )
Gdict = { 'G':G , 'ntrunc':options.lmax }
if Gdir not in glob.glob( Gdir ) :
    os.system( 'mkdir -p %s' % Gdir )
file = open( Gpath , 'wb') ; cpkl.dump( Gdict , file , -1 ) ; file.close()

days_skippedpath = Gdir + 'days_skipped.pkl'
file = open( days_skippedpath , 'wb' ) ; cpkl.dump( days_skipped , file , -1 ) ; file.close()



                            
                                                                        
