#!/usr/bin/env python
import os
import sys
import glob
import AnisotropySearch as AS
import numpy as np
import cPickle as cpkl
from optparse import OptionParser

usage = """
%prog TSDIR ORFDIR PSDDIR XDIR\n
TSDIR --- path to directory containing time-series\n
ORFDIR --- path to directory containing ORF multipole moments\n
PSDDIR --- path to directory containing PSDs\n
Xdir --- path to directory to save X in\n
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

parser.add_option( '--window' , action='store' , dest='window' , type='string' , nargs=1 ,default='None' ,
                   help='Window applied to the time-series for forming the cross-spectra' )

( options , args ) = parser.parse_args()


if len( args ) < 5 :
    parser.error( "You must specify TSDIR, PSDDIR, ORDDIR, XDIR and cIJdir! See Help: ./pklX.py -h" )

tsdir , orfdir , psddir , Xdir , cIJdir = args[ :5 ]


Xpath = Xdir + '/X.pkl'

if cIJdir not in glob.glob( cIJdir ) :
    os.system( 'mkdir -p %s' % cIJdir )

firstavailable = True ; days_skipped = []
for day in options.days :
    print "~~~ Day %d ~~~" % day
    
    tspath  = tsdir  + '/d%03d.pkl' % day
    psdpath = psddir + '/d%03d.pkl' % day
    orfpath = orfdir + '/orf_d%03d.pkl' % day

    if tspath not in glob.glob( tspath ) :
        print 'Time-series not available on day %03d' % day ; days_skipped += [ day ] ; continue
    if psdpath not in glob.glob( psdpath ) :
        print 'PSD not available on day %03d' % day ; days_skipped += [ day ] ; continue
    if orfpath not in glob.glob( orfpath ) :
        print 'ORF not found at %s' % orfpath ; days_skipped += [ day ] ; continue
#        print 'ORF not availalbe on day %03d' % day ; days_skipped += [ day ] ; continue

    print "Time-series, psd and orf all available on day %03d. " % day
    print "First available day?"  , firstavailable

    file = open( tspath , 'rb' ) ; tsdict = cpkl.load( file ) ; file.close()
    ts = AS.TimeSeries( tsdict )
    if options.scale_ts :
        ts.scale_by( options.scale_ts )

    file = open( psdpath , 'rb' ) ; psddict = cpkl.load( file ) ; file.close()
    
    orf = AS.OrfMultipleMoments( orfpath )
    
    xxx = AS.XXX( orf , psddict , ts , options.GWslope , day , cIJdir=cIJdir )
    XX = xxx.getsummand( flow=options.flow , fhigh=options.fhigh , lmax=options.lmax , window=options.window )

    if firstavailable :
        X = np.zeros( np.shape(XX.data) , complex )
        print 'Calculating normalisation factor due to coarsegraining and windowing'
        N = xxx.ts.t.data.shape[0] ; T = xxx.ts.t.Cadence1 * N ; df = xxx.fcoarse.Cadence1
        if options.window == 'None' :
            window = np.ones( N )
        elif options.window == 'hanning' :
            window = np.hanning( N ) 
        norm = ( np.sum( window**2 ) / N )**2 / ( np.sum( window**4 ) / N ) * T*df
        firstavailable = False

    dXdata = norm * XX.data
    X += dXdata

    "~~Write X_day to disk for each day here"
    map_dX = AS.xlmSkyMap( xlm = dXdata )
    if Xdir+'/../XX' not in glob.glob( Xdir+'/../XX' ) :
        os.system( 'mkdir -p %s' % ( Xdir+'/../XX' ) )
    file = open( Xdir+'/../XX/XX_d%03d.pkl' % day , 'wb' ) ; cpkl.dump( map_dX , file , -1 ) ; file.close()
    "~~"


map_X = AS.xlmSkyMap( xlm = X )

if Xdir not in glob.glob( Xdir ) :
    os.system( 'mkdir -p %s' % Xdir )
file = open( Xpath , 'wb') ; cpkl.dump( map_X , file , -1 ) ; file.close()

days_skippedpath = Xdir + 'days_skipped.pkl'
file = open( days_skippedpath , 'wb' ) ; cpkl.dump( days_skipped , file , -1 ) ; file.close()

file = open( Xdir + '/x_pklX_FIN.pkl' , 'rb' ) ; FIN = cpkl.load( file ) ; file.close()
FIN += [ 1 ] ; file = open( Xdir + '/x_pklX_FIN.pkl' , 'wb' ) ; cpkl.dump( FIN , file , -1 ) ; file.close()






