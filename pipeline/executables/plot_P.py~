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
$prog PPATH FIGPATH\n
PPATH --- path to file containing P\n
FIGPATH --- path to file to save the figure in
"""
parser = OptionParser( usage = usage )

parser.add_option( '--lmax' , action='store' , dest='lmax' , type='int' , default=15 ,
                   help='l up to which to truncate the P_{lm}s.' )

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
    parser.error( 'You must specify PPATH and FIGPATH! ' )

Ppath , figpath = args[ :2 ]



file = open( Ppath , 'rb' ) ; map_P = cpkl.load( file ) ; file.close()
Plm = AS.get_lmax_subset_from( map_P.xlm , options.lmax )
Pmap = AS.xlmSkyMap( xlm = options.mapnorm * Plm )
Pmap.xlm_to_plm() ; Pmap.xlm_to_qlm()
Pmap.create_sky( nlon=options.nlon , nlat=options.nlat )
Pmap.plm_to_P() ; Pmap.qlm_to_Q()

figdir = os.path.dirname( figpath )
if figdir not in glob.glob( figdir ) :
    os.system( 'mkdir -p %s' % figdir )
if options.plot_imag :
    Qpath = figdir + '/P_imag.png'
else :
    Qpath = None
PP.project_SkyMap( Pmap , Ppath = figpath , Qpath = Qpath )


