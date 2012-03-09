#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import optparse
import numpy as np
import AnisotropySearch as AS
import myLISAmodule as mlisar

def arbitrary_SpHreal_ML( f ) :
    g00 = 9
    SpHreal_ML = np.ones( f.shape ) * ( g00/np.sqrt(2*np.pi) )
    return SpHreal_ML

def arbitrary_SpHimag_ML( f ) :
    SpHimag_ML = np.zeros( f.shape )
    return SpHimag_ML

usage = """
%prog ORFPATH
This let's you set multipole moments of overlap-reduction function arbitrarily, and save them in the form\n
{'OrfMultipleMoments':\n
{ 'ntrunc': lmax , 'f': f , 'Antenna': IJ , 'real': SpHreal , 'imag': SpHimag } }\n
which is loadable by AnisotropySearch.OrfMultipleMoments()\n
ORFPATH --- file path to save the orfs\n 
"""
parser = optparse.OptionParser( usage = usage )

parser.add_option( '--f0' , action='store' , dest='f0' , type='string' , nargs=1 , default='.5' , help='Lowest frequency at which to evaulate ORF' )

parser.add_option( '--df' , action='store' , dest='df' , type='string' , nargs=1 , default='.1' , help='Frequency Cadence of ORF' )

parser.add_option( '--Nf' , action='store' , dest='Nf' , type='string' , nargs=1 , default='5' , help='Number of frequencies' )

parser.add_option( '--IJ' , action='store' , dest='IJ' , type='string' , nargs=1 , default='AA' , help='IJ of ORF' )

parser.add_option( '--lmax' , action='store' , dest='lmax' , type='string' , nargs=1 , default='0' , help='Maximum degree l for the SpHs' )

( options , args ) = parser.parse_args()
if len( args ) < 1 :
    parser.error( 'You must specify ORFPATH, the file path to save the tdiORF_SpHs!' )
else :
    orfpath = args[ 0 ]


f = options.f0 + options.df * np.arange( options.Nf )
indxp = AS.getMLvec( options.lmax , 'p' )

SpHreal = np.ones( ( len( indxp ) , options.Nf ) ) * arbitrary_SpHreal_ML( f )
SpHimag = np.ones( ( len( indxp ) , options.Nf ) ) * arbitrary_SpHimag_ML( f )

orfdict = {'OrfMultipleMoments': { 'ntrunc': options.lmax , 'f': f , 'Antenna': options.IJ , 'real': SpHreal , 'imag': SpHimag } }

orfdir = os.path.dirname( orfpath )
if orfdir not in glob.glob( orfdir ) :
    os.system( 'mkdir -p %s' % orfdir )
file = open( orfpath , 'wb' ) ; cpkl.dump( orfdict , file , -1 ) ; file.close()
print " Orf multipole moments saved in %s " % orfpath



