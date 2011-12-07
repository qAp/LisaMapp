#!/usr/bin/env python
import os
import sys
import cPickle as cpkl
import AnisotropySearch as AS
import glob 
import numpy as np
import optparse 

"""
This computes the expected PSDs given some P_{lm}
"""
usage = """
%prog PPATH ORFDIR PSDDIR\n
PPATH --- path to the skymap whose signal is to be simulated.\n
ORFDIR --- directory to the overlap reductions' multipole moments(effective the GW response to the skymap).\n
PSDDIR --- path to the file in which to save the expected PSDs.
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
    Ppath , orfdir , psddir = args[ :3 ]

file = open( Ppath , 'rb' ) ; skymap = cpkl.load( file ) ; file.close

if options.whichtdi == 'ordinary' :
    IIs = [ '11' , '22' , '33' ]
elif options.whichtdi == 'optimal' :
    IIs = [ 'AA' , 'EE' , 'TT' ]

if psddir not in glob.glob( psddir ) :
    os.system( 'mkdir %s' % psddir )


for day in options.days :

    print "~~~~~~~ Day %d ~~~~~~~" % day
    orfids = [ ( options.tditype , options.tdigen , II[0] , options.tditype , options.tdigen , II[1] ,
                 options.f0 , options.df , options.Nf , day ) for II in IIs ]
    orfpaths  = [ ( orfdir + 'tdiI_%s_%s_%s_tdiJ_%s_%s_%s_lmax_20_f0_%f_df_%f_Nf_%d/data_nlon_120_nlat_61/orf_d%03d.pkl'
                    % orfid ) for orfid in orfids ]
    orfs     = [ AS.OrfMultipleMoments( orfpath ) for orfpath in orfpaths ]

    PIIs = [ AS.Convolve( orf , skymap , options.GWslope ) for orf in orfs ]
    P11 , P22 , P33 = PIIs

    fdata = P11.Offset1 + P11.Cadence1 * np.arange( P11.data.shape[0] )
    f = AS.Coarsable( fdata , Offset1=P11.Offset1 , Cadence1=P11.Cadence1 )

    psddict = { 'f':f , 'AA':P11 , 'EE':P22 , 'TT':P33  }
    
    file = open( psddir + 'd%03d.pkl' % day , 'wb' ) ; cpkl.dump( psddict , file , -1 ) ; file.close()

    
