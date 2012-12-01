#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
import AnisotropySearch as AS
import PostProcess as PP

usage = """
$prog XPATH FIGPATH\n
XPATH --- path to file containing X\n
FIGPATH --- path to file to save the figure in
"""
parser = OptionParser( usage = usage )
parser.add_option( '--lmax' , action='store' , dest='lmax' , type='int' , default=15 ,
                   help='l up to which to truncate the X_{lm}s.' )
parser.add_option( '--mapnorm' , action='store' , dest='mapnorm' , type='float' , default=1. ,
                   help='Scalar factor to scale all multipole moments before plotting.' )
parser.add_option( '--nlon' , action='store' , dest='nlon' , type='int' , default=180 ,
                   help='Number of longitudes to plot in the sky.' )
parser.add_option( '--nlat' , action='store' , dest='nlat' , type='int' , default=91 ,
                   help='Number of latitudes to plot in the sky.' )
parser.add_option( '--plot_imag' , action='store_true' ,
                   help='Plot the imaginary part as well and save in the same directory as the real part.' )

( options , args ) = parser.parse_args()

if len( args ) < 2 :
    parser.error( 'You must specify XPATH and FIGPATH! ' )

Xpath , figpath = args[ :2 ]



file = open( Xpath , 'rb' ) ; map_X = cpkl.load( file ) ; file.close()
Xlm = AS.get_lmax_subset_from( map_X.xlm , options.lmax )
Xmap = AS.xlmSkyMap( xlm = options.mapnorm * Xlm )
Xmap.xlm_to_plm() ; Xmap.xlm_to_qlm()
Xmap.create_sky( nlon=options.nlon , nlat=options.nlat )
Xmap.plm_to_P() ; Xmap.qlm_to_Q()

figdir = os.path.dirname( figpath )
if figdir not in glob.glob( figdir ) :
    os.system( 'mkdir -p %s' % figdir )
if options.plot_imag :
    Qpath = figdir + '/X_imag.png'
else :
    Qpath = None
PP.project_SkyMap( Xmap , Ppath = figpath , Qpath = Qpath )


