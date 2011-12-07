#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import AnisotropySearch as AS

IJs = [ 'AE' , 'AT' , 'ET' ,  'AE_AT' , 'AE_ET' , 'AT_ET' , 'AE_AT_ET' ]
injPpath = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/skymaps/library/sphericalharmonics/Y_l0_m0_x1e-34/Y_l0_m0.pkl'
anadir = 'GW_slope_0_test1/'
m , l , lmax = 0 , 0 , 15
scale_Plm = 1.

indxpn = AS.getMLvec( lmax )

""" Load injected Plm """
file = open( injPpath ) ; P0 = cpkl.load( file ) ; file.close()


""" Load clean map of Plm """
Plms = []
for IJ in IJs :
    Ppath = '%s/%s/P/P_lmax_%d.pkl' % ( anadir , IJ , lmax )
    file = open( Ppath , 'rb' ) ; P = cpkl.load( file ) ; file.close()
    Plms += [ P.plm[ indxpn.index( (m , l) ) ] * scale_Plm ]


""" Load standard deviation of Plm """
stdPlms = []
for IJ in IJs :
    stdPpath = '%s/%s/stdP/stdP_lmax_%d.pkl' % ( anadir , IJ , lmax )
    file = open( stdPpath , 'rb' ) ; stdP = cpkl.load( file ) ; file.close()
    stdPlms += [ stdP[ indxpn.index( (m , l) ) ] * scale_Plm ]


print '%18s'*( len(IJs) + 1 ) % tuple( [ '(%d, %d, %d)' % (m , l , lmax) ] + IJs )
print ( '%18s'+'%18.8e'*len(IJs) )  % tuple( [ 'Plm' ] + list( Plms ) )
print ( '%18s'+'%18.8e'*len(IJs) )  % tuple( [ 'Plm - P0lm' ] + list( Plms - P0.plm[ indxpn.index((m , l)) ] ) )
print ( '%18s'+'%18.8e'*len( IJs ) ) % tuple( [ 'sigma_lm' ] + list( stdPlms ) )
