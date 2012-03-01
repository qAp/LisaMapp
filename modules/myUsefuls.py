import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import scipy as sp
import scipy.linalg as sp_linalg

"""
Non-LISA-specific routines
"""



def get_noise_freq_domain_1NSD( P , df , inittime , parityN , seed ) :
    """
    returns the noise time-sereis given the power spectral density
    INPUT:
    P --- power spectral density. numpy array
    df --- frequency resolution
    inittime --- initial time of the noise time-series
    parityN --- is the length of the time-series 'Odd' or 'Even'
    seed --- seed for the noise
    OUPUT:
    t --- time [s]
    n --- noise time-series
    """
    Nf = P.shape[0]
    
    if parityN == 'Odd' :
        N = 2 * Nf + 1
    elif parityN == 'Even' :
        N = 2 * ( Nf + 1 ) 
    else :
        raise InputError , "parityN must be either 'Odd' or 'Even'!"
    stime = 1 / ( N*df )
    t = inittime + stime * np.arange( N )
    
    np.random.seed( seed )
    z = ( np.random.standard_normal( P.shape ) + 1j * np.random.standard_normal( P.shape ) ) / np.sqrt( 2 )

    ntilde_fplus = np.sqrt( N / ( 2*stime ) * P ) * z 

    if N % 2 == 0 :
        ntilde = np.array( [0] + list( ntilde_fplus ) + [0] + list( np.flipud(np.conj(ntilde_fplus)) ) )
    else :
        ntilde = np.array( [0] + list( ntilde_fplus ) + list( np.flipud(np.conj(ntilde_fplus)) ) )
    n = np.real( sp.ifft( ntilde ) )
    return t , n


def get_noise_freq_domain_CovarMatrix( comatrix , df , inittime , parityN , seed='none' , N_previous_draws=0 ) :
    """
    returns the noise time-series given their covariance matrix
    INPUT:
    comatrix --- covariance matrix, Nts x Nts x Nf numpy array
          ( Nts = number of time-series. Nf number of positive and non-Nyquist frequencies )
    df --- frequency resolution
    inittime --- initial time of the noise time-series
    parityN --- is the length of the time-series 'Odd' or 'Even'
    seed --- seed for the random number generator
    N_previous_draws --- number of random number draws to discard first 
    OUPUT:
    t --- time [s]
    n --- noise time-series, Nts x N numpy array
    """
    if len( comatrix.shape ) != 3 :
        raise InputError , 'Input Covariance matrices must be a 3-D numpy array!'
    if comatrix.shape[0]  != comatrix.shape[1] :
        raise InputError , 'Covariance matrix must be square at each frequency!'

    Nts , Nf = comatrix.shape[0] , comatrix.shape[2]

    if parityN == 'Odd' :
        N = 2 * Nf + 1
    elif parityN == 'Even' :
        N = 2 * ( Nf + 1 ) 
    else :
        raise InputError , "parityN must be either 'Odd' or 'Even'!"
    stime = 1 / ( N*df )
    t = inittime + stime * np.arange( N )

    if seed == 'none' :
        print 'Not setting the seed for np.random.standard_normal()'
        pass
    elif seed == 'random' :
        np.random.seed( None )
    else :
        np.random.seed( int( seed ) )

    np.random.standard_normal( N_previous_draws ) ;

    zs = np.array( [ ( np.random.standard_normal((Nf,)) + 1j * np.random.standard_normal((Nf,)) ) / np.sqrt(2)
                     for i in range( Nts ) ] )

    ntilde_p = np.zeros( ( Nts , Nf ) , dtype=complex )
    for k in range( Nf ) :
        C = comatrix[ :,:,k ]
        if not np.allclose( C , np.conj( np.transpose( C ) ) ) :
            print "Covariance matrix NOT Hermitian! Unphysical."        
        w , V = sp_linalg.eig( C )
        for ww in w :
            if ww < 0 :
                print "Negative Eigenvalue, trying to simulate unphysical signal!"
        ntilde_p[ :,k ] =  np.conj( np.sqrt( N / (2*stime) ) * np.dot( V , np.dot( np.sqrt( np.diag( w ) ) , zs[ :,k ] ) ) )
    
    zerofill = np.zeros( ( Nts , 1 ) )
    if N % 2 == 0 :
        ntilde = np.concatenate( ( zerofill , ntilde_p , zerofill , np.conj(np.fliplr(ntilde_p)) ) , axis = 1 )
    else :
        ntilde = np.concatenate( ( zerofill , ntilde_p , np.conj(np.fliplr(ntilde_p)) ) , axis = 1 )
    n = np.real( sp.ifft( ntilde , axis = 1 ) )
    return t , n
