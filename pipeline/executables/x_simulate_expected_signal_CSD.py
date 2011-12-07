#!/usr/bin/env python
import os
import sys
import cPickle as cpkl
import AnisotropySearch as AS
import glob 
import numpy as np
import optparse

"""
This computes the expected CSDs given some P_{lm}
"""
usage = """
%prog PPATH ORFDIR CSDDIR\n
PPATH --- path to the skymap whose signal is to be simulated.\n
ORFDIR --- directory to the overlap reductions' multipole moments(effective the GW response to the skymap).\n
CSDDIR --- path to the file in which to save the expected CSDs.
"""

parser = optparse.OptionParser( usage = usage )

parser.add_option( '--GWslope' , action='store' , dest='GWslope' , type='int' , default=0 , nargs=1 ,
                   help='Gravitational spectral slope' )

parser.add_option( '--tditype' , action='store' , dest='tditype' , type='string' , default='Michelson' , nargs=1 ,
                   help='Type of TDI observable' )

parser.add_option( '--tdigen' , action='store' , dest='tdigen' , type='string' , default='G0' , nargs=1 ,
                   help='Generation of TDI observable' )

parser.add_option( '--whichtdi' , action='store' , dest='whichtdi' , type='choice' , default='optimal' , choices=('ordinary','optimal') , nargs=1 ,
                   help="Which TDI observables to simulate: 'ordinary' or 'optimal'" )

parser.add_option( '--lmax' , action='store' , dest='lmax' , type='int' , default=15 , nargs=1 , 
                   help='l up to which to convolve P and ORF' )

parser.add_option( '--f0' , action='store' , dest='f0' , type='float' , default=0.0001 , nargs=1 , 
                   help='Lower frequency of ORF' )

parser.add_option( '--df' , action='store' , dest='df' , type='float' , default=0.0001 , nargs=1 ,
                   help='Frequency resolution of ORF' )

parser.add_option( '--Nf' , action='store' , dest='Nf' , type='int' , default=4999 , nargs=1 ,
                   help='Number of frequencies of ORF' )

parser.add_option( '-d' , '--day' , action='append' , dest='days' , type='int' , nargs=1 ,
                   help='Days to simulate signals for' )

( options , args ) = parser.parse_args()

if len( args ) < 3 :
    parser.error( "You must specify PPATH, ORFDIR AND PSDDIR!  Type './x_simulate_expected_signal_PSD.py -h' for help " )
else :
    Ppath , orfdir , csddir = args[ :3 ]

file = open( Ppath , 'rb' ) ; skymap = cpkl.load( file ) ; file.close

if options.whichtdi == 'ordinary' :
    IIs = [ '12' , '13' , '23' ]
elif options.whichtdi == 'optimal' :
    IIs = [ 'AE' , 'AT' , 'ET' ]

if csddir not in glob.glob( csddir ) :
    os.system( 'mkdir %s' % csddir )


for day in options.days :

    print "~~~~~~~ Day %d ~~~~~~~" % day
    orfids = [ ( options.tditype , options.tdigen , II[0] , options.tditype , options.tdigen , II[1] ,
                 options.f0 , options.df , options.Nf , day ) for II in IIs ]
    orfpaths  = [ ( orfdir + 'tdiI_%s_%s_%s_tdiJ_%s_%s_%s_lmax_20_f0_%f_df_%f_Nf_%d/data_nlon_120_nlat_61/orf_d%03d.pkl'
                    % orfid ) for orfid in orfids ]
    orfs     = [ AS.OrfMultipleMoments( orfpath ) for orfpath in orfpaths ]

    PIJs = [ AS.Convolve( orf , skymap , options.GWslope ) for orf in orfs ]

    P12 , P13 , P23 = PIJs

    fdata = P12.Offset1 + P12.Cadence1 * np.arange( P12.data.shape[0] )
    f = AS.Coarsable( fdata , Offset1=P12.Offset1 , Cadence1=P12.Cadence1 )

    csddict = { 'f':f , 'AE':P12 , 'AT':P13 , 'ET':P23  }
    
    file = open( csddir + 'd%03d.pkl' % day , 'wb' ) ; cpkl.dump( csddict , file , -1 ) ; file.close()

    
