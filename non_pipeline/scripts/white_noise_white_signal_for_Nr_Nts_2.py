#!/usr/bin/env python
import os
import sys
import numpy as np
import cPickle as cpkl
import pylab as pl
import scipy as sp
import matplotlib.pyplot as plt
import AnisotropySearch as AS


seed = None
stime , N = .5 , 1000 ; fs = 1. / stime
var_n1 , var_n2 , Omega = 1.0 , 1.0 , 0.1
nfft = int( N / 15 ) ; noverlap = nfft / 2 ; df = 1./( nfft*stime )
numrs = range( 1 , 1000+1 )


np.random.seed( seed )

Ometers , Ometer1s , var_Ometer1s , var_Ometer1s_weak = [] , [] , [] , []
for r , numr in enumerate( numrs ) :
    n1 = np.sqrt( var_n1 ) * np.random.standard_normal( ( N, ) )
    n2 = np.sqrt( var_n2 ) * np.random.standard_normal( ( N, ) )
    h  = np.sqrt( Omega ) * np.random.standard_normal( ( N, ) )
    y1 = n1 + h ; y2 = n2 + h

    n1l = np.sqrt( var_n1 ) * np.random.standard_normal( ( N, ) )
    n2l = np.sqrt( var_n2 ) * np.random.standard_normal( ( N, ) )
    hl  = np.sqrt( Omega ) * np.random.standard_normal( ( N, ) )
    y1l = n1l + hl ; y2l = n2l + hl
    
    n1r = np.sqrt( var_n1 ) * np.random.standard_normal( ( N, ) )
    n2r = np.sqrt( var_n2 ) * np.random.standard_normal( ( N, ) )
    hr  = np.sqrt( Omega ) * np.random.standard_normal( ( N, ) )
    y1r = n1r + hr ; y2r = n2r + hr

    """ Straight cross-correlation """
    Ometer = np.mean( y1 * y2 )
    Ometer_ana , var_Ometer_ana = Omega , ( var_n1 + Omega )*( var_n2 + Omega ) / N + Omega**2/N

    """ Cross-correlation with filtering """
    Y1 = stime * np.fft.fft( y1 ) ; Y2 = stime * np.fft.fft( y2 ) ; df_c12 = 1. / ( N*stime )
    if N % 2 == 0 :
        c12 = 2 * np.conj( Y1[ : N/2+1 ] ) * Y2[ : N/2+1 ] / (N*stime) ; f_c12 = df_c12 * np.arange( N/2+1 )
    else :
        c12 = 2 * np.conj( Y1[ : (N+1)/2 ] ) * Y2[ : (N+1)/2 ] / (N*stime) ; f_c12 = df_c12 * np.arange( (N+1)/2 )
    Qc12 = AS.Coarsable( c12 , Offset1=f_c12[0] , Cadence1=df_c12 )
    Qf_c12 = AS.Coarsable( f_c12 , Offset1=f_c12[0] , Cadence1=df_c12 )

    p11l , f_p11l = pl.psd( y1l , nfft , fs , noverlap=noverlap )
    p11r , f_p11r = pl.psd( y1r , nfft , fs , noverlap=noverlap )
    p11 = ( p11l + p11r ) / 2 ; f_p11 = np.copy( f_p11l ) ; df_p11 = f_p11[1] - f_p11[0]
    Qp11 = AS.Coarsable( p11 , Offset1=f_p11[0] , Cadence1=df_p11 )
    Qf_p11 = AS.Coarsable( f_p11 , Offset1=f_p11[0] , Cadence1=df_p11 )
    p22l , f_p22l = pl.psd( y2l , nfft , fs , noverlap=noverlap )
    p22r , f_p22r = pl.psd( y2r , nfft , fs , noverlap=noverlap )
    p22 = ( p22l + p22r ) / 2 ; f_p22 = np.copy( f_p22l ) ; df_p22 = f_p22[1] - f_p22[0]
    Qp22 = AS.Coarsable( p22 , Offset1=f_p22[0] , Cadence1=df_p22 )

    Qff = AS.coarsefrequency( Qf_c12 , Qf_p11 ) ; Nf = Qff.data.shape[0]

    Qcc12 = Qc12.coarsegrain( Offset1=Qff.Offset1 , Cadence1=Qff.Cadence1 , N1=Qff.data.shape[0] )
    Qpp11 = Qp11.coarsegrain( Offset1=Qff.Offset1 , Cadence1=Qff.Cadence1 , N1=Qff.data.shape[0] )
    Qpp22 = Qp22.coarsegrain( Offset1=Qff.Offset1 , Cadence1=Qff.Cadence1 , N1=Qff.data.shape[0] )
    
    GG = (N*stime*Qff.Cadence1) * 4 * stime**2 * np.sum( 2. / (Qpp11.data * Qpp22.data) )
    XX = (N*stime*Qff.Cadence1) * 2 * stime * np.sum( 2*np.real( Qcc12.data ) / (Qpp11.data * Qpp22.data) )
    
    Ometer1 = XX / GG
    var_Ometer1 = 1. / GG * ( 1 +
                              (N*stime*Qff.Cadence1) * (2*stime)**4 * Omega**2 * np.sum( 2. / (Qpp11.data**2 * Qpp22.data**2) ) / GG )
    var_Ometer1_weak = 1. / GG 
    Ometers += [ Ometer ] ; Ometer1s += [ Ometer1 ]  ; var_Ometer1s += [ var_Ometer1 ] ; var_Ometer1s_weak += [ var_Ometer1_weak ]

Ometers = np.array( Ometers ) ; Ometer1s = np.array( Ometer1s ) ; var_Ometer1s = np.array( var_Ometer1s )
Ometer1_opt = np.sum( Ometer1s / var_Ometer1s ) / np.sum( 1. / var_Ometer1s )

var_Ometer1_ana_0 = ( var_n1 + Omega )*( var_n2 + Omega ) / N
var_Ometer1_ana_h = Omega**2 / N
var_Ometer1_ana = var_Ometer1_ana_0 + var_Ometer1_ana_h



""" Pack results """


results = { 'Ometers':Ometers , 'Ometer1s':Ometer1s , 'var_Ometer1s':var_Ometer1s , 'var_Ometer1s_weak':var_Ometer1s_weak , 'Ometer1_opt':Ometer1_opt ,
            'Omega':Omega , 'var_Ometer1_ana':var_Ometer1_ana ,
            'numrs':np.array( numrs ) }

file = open( 'results.pkl' , 'wb' ) ; cpkl.dump( results , file , -1 ) ; file.close()


#print Ometer1
#print Ometer1 - Omega
#print np.sqrt( var_Ometer1 )
#print np.sqrt( var_Ometer1_ana )
