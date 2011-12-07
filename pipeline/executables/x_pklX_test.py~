#!/usr/bin/env python
import os
import sys
import glob
import AnisotropySearch as AS
import numpy as np
import cPickle as cpkl
from optparse import OptionParser

usage = """
%prog TSDIR CSDDIR ORFDIR PSDDIR XDIR\n
TSDIR --- path to directory containing the time-series\n
CSDDIR --- path to directory containing cross-spectra data or CSDs\n
ORFDIR --- path to directory containing ORF multipole moments\n
PSDDIR --- path to directory containing PSDs\n
Xdir --- path to directory to save X in
cIJdir --- path to directory to save cIJs in
"""
parser = OptionParser( usage = usage )

parser.add_option( '--scale_ts' , action='store' , dest='scale_ts' , type='float' , default=1 , nargs=1 ,
                   help='Scale the time-series beforehand.' )

parser.add_option( '-d' , '--day' , action='append' , dest='days' , type='int' , nargs=1 ,
                   help='Days to sum over for maximum likelihood estimator of X' )

parser.add_option( '--GWslope' , action='store' , dest='GWslope' , type='int' , default=0 , nargs=1 ,
                   help='Gravitational spectral slope' )

parser.add_option( '--flow' , action='store' , dest='flow' , type='float' , default=1e-4 , nargs=1 ,
                   help='Lower limit of freqeuncy range for maximum likelihood estimator of X' )

parser.add_option( '--fhigh' , action='store' , dest='fhigh' , type='float' , default=1e-1 , nargs=1 ,
                   help='Higher limit of frequency range for maximum likelihood estimator of X' )

parser.add_option( '--lmax' , action='store' , dest='lmax' , type='int' , default=20 , nargs=1 ,
                   help='Maximum l to take into account in estimating X' )

( options , args ) = parser.parse_args()


if len( args ) < 5 :
    parser.error( "You must specify TSDIR, PSDDIR, ORDDIR and XDIR! See Help: ./pklX.py -h" )

tsdir , csddir , orfdir , psddir , Xdir , cIJdir = args[ :6 ]

Xpath = Xdir + '/X.pkl'

if cIJdir not in glob.glob( cIJdir ) :
    os.system( 'mkdir %s' % cIJdir )


firstavailable = True ; days_skipped = []
for day in options.days :
    print "~~~ Day %d ~~~" % day

    tspath = tsdir + '/d%03d.pkl' % day
    csdpath = csddir + '/d%03d.pkl' % day
    psdpath = psddir + '/d%03d.pkl' % day
    orfpath = orfdir + '/orf_d%03d.pkl' % day

    if tspath not in glob.glob( tspath ) :
        print 'Time-series data not availabe on day %03d' % day ; days_skipped += [ day ] ; continue
    if csdpath not in glob.glob( csdpath ) :
        print 'Cross-spectra data or CSDs not available on day %03d' % day ; days_skipped += [ day ] ; continue
    if psdpath not in glob.glob( psdpath ) :
        print 'PSD not available on day %03d' % day ; days_skipped += [ day ] ; continue
    if orfpath not in glob.glob( orfpath ) :
        print 'ORF not availalbe on day %03d' % day ; days_skipped += [ day ] ; continue

    print "ts, csd, psd and orf all available on day %03d. " % day
    print "First available day?"  , firstavailable

    file = open( csdpath , 'rb' ) ; csddict = cpkl.load( file ) ; file.close()

    file = open( psdpath , 'rb' ) ; psddict = cpkl.load( file ) ; file.close()
    
    orf = AS.OrfMultipleMoments( orfpath )
    
    xxx = AS.XXX_test( orf , psddict , csddict , options.GWslope , day , cIJdir=cIJdir )
    XX = xxx.getsummand( flow=options.flow , fhigh=options.fhigh , lmax=options.lmax )

    if firstavailable :
        X = np.zeros( np.shape(XX.data) , complex )
        firstavailable = False

    X += XX.data

file = open( tspath , 'rb' ) ; tsdict = cpkl.load( file ) ; file.close()
'normalise for coarsegraining and windowing with a Hanning window for both time-series'
ts = AS.TimeSeries( tsdict )
N = ts.t.data.shape[0] ; T = ts.t.Cadence1 * N ; df = xxx.fcoarse.Cadence1
window = np.hanning( N )
norm = ( np.sum( window**2 ) / N ) / ( np.sum( window**4 ) / N ) * T*df
map_X = AS.xlmSkyMap( xlm = norm * X )

if Xdir not in glob.glob( Xdir ) :
    os.system( 'mkdir -p %s' % Xdir )
file = open( Xpath , 'wb') ; cpkl.dump( map_X , file , -1 ) ; file.close()

days_skippedpath = Xdir + '/days_skipped.pkl'
file = open( days_skippedpath , 'wb' ) ; cpkl.dump( days_skipped , file , -1 ) ; file.close()




