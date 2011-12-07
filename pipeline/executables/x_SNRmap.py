#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import  optparse
import AnisotropySearch as AS

usage = """
%prog GPATH XPATH SNRMAPPATH\n
GPATH   --- path to file containing Fisher matrix G\n
XPATH   --- path to file containing dirty map X\n
SNRMAPPATH --- path to file in which to save sigmamap\n
"""
parser = optparse.OptionParser( usage = usage )
parser.add_option( '--regMethod' , action='store' , dest='regMethod' , type='int' , default=1 ,
                   help='The method to regularise the Fisher matrix with.' )

parser.add_option( '--regCutoff' , action='store' , dest='regCutoff' , type='string' , default='9.5e-4' ,
                   help='The cut-off for the regularisation of the Fisher matrix.' )

parser.add_option( '--nlon' , action='store' , dest='nlon' , type='int' , default=180 ,
                   help='Number of longitudes to plot in the sky.' )

parser.add_option( '--nlat' , action='store' , dest='nlat' , type='int' , default=91 ,
                   help='Number of latitudes to plot in the sky.' )

parser.add_option( '--lmax' , action='store' , dest='lmax' , type='int' , default=15 ,
                   help='l up to which to truncate the X and G.' )

parser.add_option( '--mapnorm' , action='store' , dest='mapnorm' , type='string' , default='1.' ,
                   help='Scalar factor to scale the map plotting.' )

( options , args ) = parser.parse_args()

if len( args ) < 3 :
    parser.error( 'You must specify at least GPATH, XPATH and SNRMAPPATH! Type ./x_SNRmap.py -h' )

Gpath , Xpath , SNRmappath = args[ :3 ]

fish = AS.FisherMatrix( Gpath , lmax=options.lmax )
if options.regMethod == 0 :
    print 'Calculating unregularised inverse of Fisher matrix...'
    fishinv = fish.invert()
    print 'done'
else :
    print "Calculating regularised inverse of Fisher matrix..."
    fish.regularise( regMethod=options.regMethod , regCutoff=float( options.regCutoff ) )
    fishinv = fish.reginvert()

file = open( Xpath , 'rb' ) ; map_X = cpkl.load( file ) ; file.close()
Xlm = AS.get_lmax_subset_from( map_X.xlm , options.lmax )

Pdata = np.dot( fishinv , Xlm )
Pmap = AS.xlmSkyMap( xlm = Pdata )
Pmap.xlm_to_plm() ; Pmap.xlm_to_qlm()
Pmap.create_sky( nlon = options.nlon , nlat = options.nlat )
Pmap.plm_to_P() ; Pmap.qlm_to_Q() ; Pmap.PQ_to_X()

covarm = np.dot( fishinv , np.dot( fish.fish , fishinv ) )
map_sigma , lats , lons = AS.getSigmaMap( covarm , options.nlat , options.nlon )
mapdata = np.reshape( map_sigma , ( options.nlat , options.nlon ) )
sigmamap = AS.XSkyMap( X = mapdata ) ; sigmamap.X_to_P() ; sigmamap.X_to_Q()
sigmamap.ntrunc_equals( options.lmax ) ; sigmamap.P_to_plm() ; sigmamap.Q_to_qlm() ; sigmamap.plmqlm_to_xlm()

SNRdata = Pmap.P / sigmamap.P
SNRmap = AS.XSkyMap( X = float( options.mapnorm ) * SNRdata ) ; SNRmap.X_to_P() ; SNRmap.X_to_Q()
SNRmap.ntrunc_equals( options.lmax ) ; SNRmap.P_to_plm() ; SNRmap.Q_to_qlm() , SNRmap.plmqlm_to_xlm()

SNRmapdir = os.path.dirname( SNRmappath )
if SNRmapdir not in glob.glob( SNRmapdir ) :
    os.system( 'mkdir %s' % SNRmapdir )
file = open( SNRmappath , 'wb' ) ; cpkl.dump( SNRmap , file , -1 ) ; file.close()



