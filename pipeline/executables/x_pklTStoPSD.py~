#!/usr/bin/env python
import os
import glob
import sys
import cPickle as cpkl
import AnisotropySearch as AS
import pylab as pl
import optparse



parser = optparse.OptionParser( 'Usage: ./x_pklTStoPSD.py TSDIR PSDDIR' )
parser.add_option( '--scale_ts' , action='store' , dest='scale_ts' , type='float' , default=1 , nargs=1 ,
                   help='Scale the time-series beforehand.' )
parser.add_option( '-d' , '--day' , action='append' , dest='days' , type='int' , nargs = 1 ,
                   help='Days on which to estimate PSD.' )
parser.add_option( '--segduration' , action='store' , dest='segduration' , type='float' , default=2*60**2 ,
                   help='Duration of averaging segments in Pwelch PSD estimation[s]' )

( options , args ) = parser.parse_args()

if len( args ) < 2 :
    parser.error( "You must specify TSDIR, directory containing time-series, and PSDDIR, the directory to save the PSDs in!" )

tsdir  = args[0]
scale_ts = options.scale_ts
days = options.days
segduration = options.segduration
psddir = args[1]



if psddir not in glob.glob( psddir ) :
    os.system( 'mkdir %s' % psddir )

for day in options.days :
    print "---- DAY %d ----" % day
    tspath = tsdir + '/d%03d.pkl' % day
    if tspath not in glob.glob( tspath ) :
        print 'Time-series NOT found at %s' % tspath
        continue

    psdpath = psddir + 'd%03d.pkl' % day
    if psdpath in glob.glob( psdpath ) :
        print 'psd/%s exists. Nothing to do.' % psdpath
        continue

    print "Time-series available on day %03d, calculating its PSDs... " % day ,

    file = open( tspath , 'rb' ) ; tsdict = cpkl.load( file ) ; file.close()

    ts = AS.TimeSeries( tsdict )
    if options.scale_ts :
        ts.scale_by( options.scale_ts )

    t , n1 , n2 , n3 = ts.t.data , ts.A.data , ts.E.data , ts.T.data
    stime = ts.t.Cadence1 ; inittime = ts.t.Offset1 ; N = t.shape[0]

    nfft = int( options.segduration / stime ) ; noverlap = nfft / 2 ; fs = 1. / stime

    P11data , fdata = pl.psd( n1 , nfft , fs , noverlap = noverlap )
    P22data , fdata = pl.psd( n2 , nfft , fs , noverlap = noverlap )
    P33data , fdata = pl.psd( n3 , nfft , fs , noverlap = noverlap )

    f0 = fdata[0] ; df = fdata[1] - fdata[0] ; fscale = { 'Offset1':f0 , 'Cadence1':df }
    f = AS.Coarsable( fdata , **fscale )
    P11 = AS.Coarsable( P11data , **fscale )
    P22 = AS.Coarsable( P22data , **fscale )
    P33 = AS.Coarsable( P33data , **fscale )

    psddict = { 'f': f , 'AA': P11 , 'EE': P22 , 'TT': P33 }



    file = open( psdpath , 'wb' )
    cpkl.dump( psddict , file , -1 )
    file.close()
    
    print "done"



