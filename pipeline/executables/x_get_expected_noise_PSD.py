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
This computes the expected noise PSDs 
"""
usage = """
%prog PSDDIR\n
PSDDIR --- path to the file in which to save the expected PSDs.
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
    parser.error( "You must specify PSDDIR!  Type './x_simulate_expected_signal_PSD.py -h' for help " )
else :
    psddir = args[ 0 ]

if options.whichtdi == 'ordinary' :
    IIs = [ '11' , '22' , '33' ]
elif options.whichtdi == 'optimal' :
    IIs = [ 'AA' , 'EE' , 'TT' ]

if psddir not in glob.glob( psddir ) :
    os.system( 'mkdir %s' % psddir )

fdata = options.f0 + options.df * np.arange( options.Nf )
fscale = { 'Offset1':options.f0 , 'Cadence1':options.df }
f = AS.Coarsable( fdata , **fscale )

print "Calculating expected auto NSD..."
PIIdatas = [ mlisar.get_tdiNSD( options.tditype , options.tdigen , II[0] , II[1] , fdata ) for II in IIs ]
PIIs = [ AS.Coarsable( PIIdata , **fscale ) for PIIdata in PIIdatas ]
P11 , P22 , P33 = PIIs
print 'done'

psddict = { 'f':f , 'AA':P11 , 'EE':P22 , 'TT':P33  }
    
file = open( psddir + 'd001.pkl' , 'wb' ) ; cpkl.dump( psddict , file , -1 ) ; file.close()

    
