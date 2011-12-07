#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import matplotlib.pyplot as plt
import myLISAmodule as mlisar
import AnisotropySearch as AS





t0s = [ 100 ]
m , l = 3 , 5
tdiORF_SpH_dir = 'data_nlon_220_nlat_111/'
fig_dir = 'figures_f_vs_tdiORF_SpHs_nlon_220_nlat_111/' 


indxpn = AS.getMLvec( l , 'pn' ) ; k = indxpn.index( (m,l) )

for t0 in t0s :

    orf = AS.OrfMultipleMoments( tdiORF_SpH_dir + 'orf_t0_%.1f.pkl' % t0 )
    glm = orf.getMultipleMoments( 'pn' , l )
    f = glm.Offset1 + glm.Cadence1 * np.arange( glm.data.shape[1] )
    glmdata = glm.data[ k , : ]

    fig_dir = fig_dir + 't0_%.1f/' % t0
    if fig_dir not in glob.glob( fig_dir ) :
        os.system( 'mkdir -p %s' % fig_dir )

    fig = plt.figure()
    fig.suptitle( 'tdiORF SpH: l = %d , m = %d , t0 = %.1f' % ( l , m , t0 ) )
    ax = fig.add_subplot( 311 )
    ax.plot( f , np.real( glmdata ) , f , np.imag( glmdata ) )
    ax.legend( ( 'real' , 'imag' ) )
    ax.set_ylabel( 'glm' ) ; ax.set_xlabel( 'frequency [Hz]' )
    fig.savefig( fig_dir + 'l_%03d_m_%+04d.png' % ( l , m ) )
    
