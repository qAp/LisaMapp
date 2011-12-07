#!/usr/bin/env python
import os
import sys
import cPickle as cpkl
import AnisotropySearch as AS
import glob 
import numpy as np
import optparse 
import myLISAmodule as mlisar

"""
This computes the expected noise CSDs 
"""
usage = """
%prog CSDDIR\n
CSDDIR --- path to the file in which to save the expected CSDs.
"""
parser = optparse.OptionParser( usage = usage )


parser.add_option( '--tditype' , action='store' , dest='tditype' , type='string' , default='Michelson' , nargs=1 ,
                   help='Type of TDI observable' )

parser.add_option( '--tdigen' , action='store' , dest='tdigen' , type='string' , default='G0' , nargs=1 ,
                   help='Generation of TDI observable' )

parser.add_option( '--whichtdi' , action='store' , dest='whichtdi' , type='choice' , default='optimal' , choices=('ordinary','optimal') , nargs=1 ,
                   help="Which TDI observables to simulate: 'ordinary' or 'optimal'" )

parser.add_option( '--f0' , action='store' , dest='f0' , type='float' , default=0.0001 , nargs=1 , 
                   help='Lower frequency of ORF' )

parser.add_option( '--df' , action='store' , dest='df' , type='float' , default=0.0001 , nargs=1 ,
                   help='Frequency resolution of ORF' )

parser.add_option( '--Nf' , action='store' , dest='Nf' , type='int' , default=4999 , nargs=1 ,
                   help='Number of frequencies of ORF' )

( options , args ) = parser.parse_args()

if len( args ) < 1 :
    parser.error( "You must specify CSDDIR!  Type './x_simulate_expected_signal_CSD.py -h' for help " )
else :
    csddir = args[ 0 ]

if options.whichtdi == 'ordinary' :
    IJs = [ '12' , '13' , '23' ]
elif options.whichtdi == 'optimal' :
    IJs = [ 'AE' , 'AT' , 'ET' ]

if csddir not in glob.glob( csddir ) :
    os.system( 'mkdir %s' % csddir )

fdata = options.f0 + options.df * np.arange( options.Nf )
fscale = { 'Offset1':options.f0 , 'Cadence1':options.df }
f = AS.Coarsable( fdata , **fscale )

print "Calculating expected cross NSD..."
PIJdatas = [ mlisar.get_tdiNSD( options.tditype , options.tdigen , IJ[0] , IJ[1] , fdata ) for IJ in IJs ]
PIJs = [ AS.Coarsable( PIJdata , **fscale ) for PIJdata in PIJdatas ]
P12 , P13 , P23 = PIJs
print 'done'

psddict = { 'f':f , 'AE':P12 , 'AT':P13 , 'ET':P23  }
    
file = open( csddir + 'd001.pkl' , 'wb' ) ; cpkl.dump( psddict , file , -1 ) ; file.close()

    
