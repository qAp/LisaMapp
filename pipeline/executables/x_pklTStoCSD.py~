#!/usr/bin/env python
import os
import glob
import sys
import cPickle as cpkl
import AnisotropySearch as AS
import pylab as pl
from optparse import OptionParser



parser = OptionParser( "usage: ./x_pklTS_CSD.py TSDIR CSDDIR" )
parser.add_option( '--scale_ts' , action='store' , dest='scale_ts' , type='float' , default=1 , nargs=1 ,
                   help='Scale the time-series beforehand.' )
parser.add_option( '-d' , '--day' , action='append' , dest='days' , type='int' , nargs = 1 ,
                   help='Days on which to estimate CSD.' )
parser.add_option( '--segduration' , action='store' , dest='segduration' , type='float' , default=2*60**2 ,
                   help='Duration of averaging segments in Pwelch CSD estimation[s]' )

( options , args ) = parser.parse_args()

if len( args ) < 2 :
    parser.error( "You must specify TSDIR, directory containing time-series, and CSDDIR, the directory to save the CSDs in!" )

tsdir , csddir  = args


if csddir not in glob.glob( csddir ) :
    os.system( 'mkdir %s' % csddir )

for day in options.days :
    print "---- DAY %d ----" % day

    tspath = tsdir + 'd%03d.pkl' % day
    if tspath not in glob.glob( tspath ) :
        print 'Time-series not found on day %d. Nothing to do.' % day
        continue

    csdpath = csddir + 'd%03d.pkl' % day
    if csdpath in glob.glob( csdpath ) :
        print '%s exists. Nothing to do.' % csdpath
        continue
    
    print "Time-series available on day %03d, calculating its CSDs... " % day ,

    file = open( tspath , 'rb' ) ; tsdict = cpkl.load( file ) ; file.close()

    ts = AS.TimeSeries( tsdict )
    if options.scale_ts :
        ts.scale_by( options.scale_ts )

    t , n1 , n2 , n3 = ts.t.data , ts.A.data , ts.E.data , ts.T.data
    stime = ts.t.Cadence1 ; inittime = ts.t.Offset1 ; N = t.shape[0]

    nfft = int( options.segduration / stime ) ; noverlap = nfft / 2 ; fs = 1. / stime

    P12data , fdata = pl.csd( n1 , n2 , nfft , fs , noverlap = noverlap )
    P13data , fdata = pl.csd( n1 , n3 , nfft , fs , noverlap = noverlap )
    P23data , fdata = pl.csd( n2 , n3 , nfft , fs , noverlap = noverlap )

    f0 = fdata[0] ; df = fdata[1] - fdata[0] ; fscale = { 'Offset1':f0 , 'Cadence1':df }
    f = AS.Coarsable( fdata , **fscale )
    P12 = AS.Coarsable( P12data , **fscale )
    P13 = AS.Coarsable( P13data , **fscale )
    P23 = AS.Coarsable( P23data , **fscale )

    csddict = { 'f': f , 'AE': P12 , 'AT': P13 , 'ET': P23 }

    file = open( csdpath , 'wb' )
    cpkl.dump( csddict , file , -1 )
    file.close()
    
    print "done"



