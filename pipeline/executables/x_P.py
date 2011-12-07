#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np
from optparse import OptionParser
import AnisotropySearch as AS

usage = """
$prog GPATH XPATH PPATH\n
GPATH --- path to file for Fisher matrix G\n
XPATH --- path to file for dirty map X\n
PPATH --- path to file to save clean map P\n
"""
parser = OptionParser( usage = usage )

parser.add_option( '--regMethod' , action='store' , dest='regMethod' , type='int' , default=1 ,
                   help='The method to regularise the Fisher matrix with.' )

parser.add_option( '--regCutoff' , action='store' , dest='regCutoff' , type='string' , default='9.5e-4' ,
                   help='The cut-off for the regularisation of the Fisher matrix.' )

parser.add_option( '--lmax' , action='store' , dest='lmax' , type='int' , default=15 ,
                   help='l up to which to truncate the X and G.' )

parser.add_option( '--mapnorm' , action='store' , dest='mapnorm' , type='float' , default=1. ,
                   help='Scalar factor to scale all multipole moments before plotting.' )

parser.add_option( '--nlon' , action='store' , dest='nlon' , type='int' , default=180 ,
                   help='Number of longitudes to plot in the sky.' )

parser.add_option( '--nlat' , action='store' , dest='nlat' , type='int' , default=91 ,
                   help='Number of latitudes to plot in the sky.' )

parser.add_option( '--N_keptSV' , action='store_true' , dest='N_keptSV' ,
                   help='Show the number of survived eigenvalues.' )

( options , args ) = parser.parse_args()

if len( args ) < 3 :
    parser.erro( 'You must specify GPATH, XPATH AND PPATH! (See help: x_P.py -h)' )

Gpath , Xpath , Ppath = args[ :3 ]


fish = AS.FisherMatrix( Gpath , lmax=options.lmax )
if options.regMethod == 0 :
    print 'Calculating unregularised inverse of Fisher matrix...'
    fishinv = fish.invert()
else :
    print "Calculating regularised inverse of Fisher matrix..."
    fish.regularise( regMethod=options.regMethod , regCutoff=float( options.regCutoff ) )
    fishinv = fish.reginvert()

if options.N_keptSV : 
    print 'Number of eigenvalues kept = ' , fish.N_keptSV

file = open( Xpath , 'rb' ) ; map_X = cpkl.load( file ) ; file.close()
Xlm = AS.get_lmax_subset_from( map_X.xlm , options.lmax )

Pdata = np.dot( fishinv , Xlm )

Pmap = AS.xlmSkyMap( xlm = options.mapnorm * Pdata )
Pmap.xlm_to_plm() ; Pmap.xlm_to_qlm()
Pmap.create_sky( nlon = options.nlon , nlat = options.nlat )
Pmap.plm_to_P() ; Pmap.qlm_to_Q() ; Pmap.PQ_to_X()

Pdir = os.path.dirname( Ppath )
if Pdir not in glob.glob( Pdir ) :
    os.system( 'mkdir -p %s' % Pdir )
file = open( Ppath , 'wb' ) ; cpkl.dump( Pmap , file , -1 ) ; file.close()




    
    

