#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np
from optparse import OptionParser
import AnisotropySearch as AS

usage = """
%prog GPATH SIGMAPPATH\n
GPATH   --- path to file for Fisher matrix G\n
SIGMAPPATH --- path to file in which to save sigmamap
"""
parser = OptionParser( usage = usage )

parser.add_option( '--regMethod' , action='store' , dest='regMethod' , type='int' , default=1 ,
                   help='The method to regularise the Fisher matrix with.' )

parser.add_option( '--regCutoff' , action='store' , dest='regCutoff' , type='string' , default='9.5e-4' ,
                   help='The cut-off for the regularisation of the Fisher matrix.' )

parser.add_option( '--nlon' , action='store' , dest='nlon' , type='int' , default=180 ,
                   help='Number of longitudes in the sky.' )

parser.add_option( '--nlat' , action='store' , dest='nlat' , type='int' , default=91 ,
                   help='Number of latitudes in the sky.' )

parser.add_option( '--lmax' , action='store' , dest='lmax' , type='int' , default=15 ,
                   help='l up to which to truncate the X and G.' )

parser.add_option( '--mapnorm' , action='store' , dest='mapnorm' , type='string' , default='1.' ,
                   help='Scalar factor to scale all multipole moments before plotting.' )

( options , args ) = parser.parse_args()

if len( args ) < 2 :
    parser.error( 'You must specify GPATH and SIGMAPPATH! See Help: ./x_sigmamap.py -h' )

Gpath , sigmappath = args[ :2 ]

fish = AS.FisherMatrix( Gpath , lmax = options.lmax )
fish.regularise( regMethod = options.regMethod , regCutoff = float( options.regCutoff ) )
regfishinv = fish.reginvert()
covarm = np.dot( regfishinv , np.dot( fish.fish , regfishinv ) )

map_sigma , lats , lons = AS.getSigmaMap( covarm , options.nlat , options.nlon )
mapdata = np.reshape( map_sigma , ( options.nlat , options.nlon ) )
sigmamap = AS.XSkyMap( X = float( options.mapnorm ) * mapdata ) ; sigmamap.X_to_P() ; sigmamap.X_to_Q()
sigmamap.ntrunc_equals( options.lmax ) ; sigmamap.P_to_plm() ; sigmamap.Q_to_qlm() ; sigmamap.plmqlm_to_xlm()

sigmapdir = os.path.dirname( sigmappath )
if sigmapdir not in glob.glob( sigmapdir ) :
    os.system( 'mkdir %s' % sigmapdir )
file = open( sigmappath , 'wb' ) ; cpkl.dump( sigmamap , file , -1 ) ; file.close()



