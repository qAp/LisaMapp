#!/usr/bin/env python
import os
import sys
import glob
import AnisotropySearch as AS
import numpy as np
import cPickle as cpkl
import optparse

usage="""
%prog TSDIR CSDDIR ORFDIR PSDDIR SPATH\n
TSDIR --- path to directory containing the time-series (for normalisation purpose only)\n
CSDDIR --- path to directory containing cross-spectra data or CSDs\n
ORFDIR --- path to directory containing ORF multipole moments\n
PSDDIR --- path to directory containing PSDs\n
SPATH --- path to file to save the bias matrix in
"""
parser = optparse.OptionParser( usage = usage )

parser.add_option( '-d' , '--day' , action='append' , dest='days' , type='int' , nargs=1 ,
                   help='Days to sum over' )

parser.add_option( '--GWslope' , action='store' , dest='GWslope' , type='int' , nargs=1 , default=0 ,
                   help='GW spectral slope' )

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
    parser.error( "You must specify TSDIR, CSDDIR , ORFDIR, PSDDIR and SPATH! Type './s_S.py -h'" )
else :
    tsdir , csddir , orfdir , psddir , Spath = args[ :5 ]

Sdir = os.path.dirname( Spath )

firstavailable , days_skipped = True , []
for day in options.days :
    print '~~~~~~~~~~ Day %d ~~~~~' % day

    csdpath = csddir + '/d%03d.pkl' % day
    orfpath = orfdir + '/orf_d%03d.pkl' % day
    psdpath = psddir + '/d%03d.pkl' % day

    if csdpath not in glob.glob( csdpath ) :
        print 'CSD not found at %s' % csdpath ; continue
    if orfpath not in glob.glob( orfpath ) :
        print 'Orf not found at %s' % orfpath ; continue
    if psdpath not in glob.glob( psdpath ) :
        print 'PSD not found at %s' % psdpath ; continue

    orf = AS.OrfMultipleMoments( orfpath )
    file = open( csdpath , 'rb' ) ; csddict = cpkl.load( file ) ; file.close()
    file = open( psdpath , 'rb' ) ; psddict = cpkl.load( file ) ; file.close()
    
    SSdata = AS.get_covariance_bias_matrix_for_the_day( orf , psddict , csddict , options.GWslope , day , options.flow , options.fhigh , options.lmax )
    
    if firstavailable :
        Sdata = np.zeros( SSdata.shape , dtype = SSdata.dtype )
        print 'Calculating normalisation factor due to coarsegraining and windowing'
        tspath = tsdir + '/d%03d.pkl' % day
        file = open( tspath , 'rb' ) ; tsdict = cpkl.load( file ) ; file.close() ; ts = AS.TimeSeries( tsdict )
        N = ts.t.data.shape[0] ; T = ts.t.Cadence1 * N 
        fcoarse = AS.coarsefrequency( orf.f , psddict['f'] , csddict['f'] ) ; df = fcoarse.Cadence1
        if options.window == 'None' :
            window = np.ones( N )
        elif options.window == 'hanning' :
            window = np.hanning( N )
        norm = ( np.sum( window**2 ) / N )**2 / ( np.sum( window**4 ) / N ) * T*df        
        firstavailable = False

    SSdata *= norm
    Sdata += SSdata

    "~~Write S_day to disk for each day here"
    SS = AS.Coarsable( SSdata ) ; SSdict = { 'S':SS , 'ntrunc':options.lmax }
    if Sdir+'/../SS' not in glob.glob( Sdir+'/../SS' ) :
        os.system( 'mkdir -p %s' % ( Sdir+'/../SS' ) )
    file = open( Sdir+'/../SS/SS_d%03d.pkl' % day , 'wb' ) ; cpkl.dump( SSdict , file , -1 ) ; file.close()
    "~~"

S = AS.Coarsable( Sdata )
Sdict = { 'S':S , 'ntrunc':options.lmax }


if Sdir not in glob.glob( Sdir ) :
    os.system( 'mkdir -p %s' % Sdir )
file = open( Spath , 'wb' ) ; cpkl.dump( Sdict , file , -1 ) ; file.close()

days_skippedpath = Sdir + '/days_skipped.pkl'
file = open( days_skippedpath , 'wb' ) ; cpkl.dump( days_skipped , file , -1 ) ; file.close()
print 'done'


