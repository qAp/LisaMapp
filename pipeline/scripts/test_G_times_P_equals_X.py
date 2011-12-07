#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import AnisotropySearch as AS
import PostProcess as PP


GWslopes = [ 0 ]
IJs = [ 'AE' , 'AT'  , 'ET' , 'AE_AT' , 'AE_ET' , 'AT_ET' , 'AE_AT_ET' ]
lmax , nlon , nlat = 15 , 180 , 91
Ppath = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/skymaps/library/sphericalharmonics/Y_l0_m0_x1e-34/Y_l0_m0.pkl'
GPnorm = 1e-40

file = open( Ppath , 'rb' ) ; P0 = cpkl.load( file ) ; file.close()
P0lm = AS.get_lmax_subset_from( np.copy( P0.xlm ) , lmax )

workdir = os.getcwd() + '/'

for slope in GWslopes :
    for IJ in IJs :
        print 'GWslope = %d , IJ = %s' % ( slope , IJ )
        Gpath = workdir + 'GW_slope_%d/%s/G/G.pkl' % ( slope , IJ )
        fish = AS.FisherMatrix( Gpath , lmax = lmax )
        Xlm = np.dot( fish.fish , P0lm )
        X0path = workdir + 'GW_slope_%d/%s/X/X.pkl' % ( slope , IJ )
        file = open( X0path ) ; X0 = cpkl.load( file ) ; file.close()
        X0lm = AS.get_lmax_subset_from( np.copy( X0.xlm ) , lmax )
        print 'The obtained X is the same as the product of the obtained Fisher mastrix and the injected P.' , np.allclose( Xlm , X0lm )
        print 'X0lm.shape , Xlm.shape = '  , X0lm.shape , ',' , Xlm.shape 
        print 'X0lm[ :10 ]'
        print X0lm[ :10 ]
        print 'Xlm[ :10 ]'
        print Xlm[ :10 ]
        X = AS.xlmSkyMap( xlm = GPnorm * Xlm )
        X.xlm_to_plm() ; X.xlm_to_qlm()
        X.create_sky( nlon = nlon , nlat = nlat )
        X.plm_to_P() ; X.qlm_to_Q() ; X.PQ_to_X()
        figpath = workdir + 'GW_slope_%d/%s/figures_X/GP_real.png' % ( slope , IJ )
        figdir = os.path.dirname( figpath )
        if figdir not in glob.glob( figdir ) :
            os.system( 'mkdir -p %s' % figdir )

        PP.project_SkyMap( X , Ppath = figpath )

