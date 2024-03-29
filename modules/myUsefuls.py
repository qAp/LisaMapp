import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import scipy as sp
import scipy.linalg as sp_linalg
import spharm
import AnisotropySearch as AS

"""
Non-LISA-specific routines
"""

def get_freqs_from_duration_and_stime( stime , duration ) :
    """
    Returns the positive frequencies (without DC and Nyquist) of the
    Fourier transform of a time-series with some given duration and stime.
    INPUT:
    stime --- sampling time [s]
    duration --- desired duration of time-series [s]
    OUTPUT:
    freqdict --- dictionary with the following fields:
    N --- number of samples in the time-series
    dt --- sampling time
    T --- duration of time-series
    parityN --- Whether N is an even or odd number
    Nf --- number of positive frequencies (without DC and Nyquist) of F.T
    df --- frequency resolution of the Fourier transform (F.T)
    f --- numpy array containing these frequencies
    """
    N = np.round( duration / stime )
    df = 1 / ( N * stime )
    if N % 2 == 0 :
        Nf = int( N/2 - 1 ) ; parityN = 'Even'
    else :
        Nf = int( (N-1) / 2 ) ; parityN = 'Odd'
    f = df * np.arange( 1 , Nf + 1 )
    freqdict = { 'N':N , 'dt':stime , 'T':N*stime , 'parityN':parityN ,
                 'Nf':Nf , 'df':df , 'f':f }
    return freqdict



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
    print N_previous_draws
    np.random.standard_normal( N_previous_draws ) ;

    zs = np.array( [ ( np.random.standard_normal((Nf,)) + 1j * np.random.standard_normal((Nf,)) ) / np.sqrt(2)
                     for i in range( Nts ) ] )

    ntilde_p = np.zeros( ( Nts , Nf ) , dtype=complex )
    for k in range( Nf ) :
        C = comatrix[ :,:,k ]
        if not np.allclose( C , np.conj( np.transpose( C ) ) ) :
            print "Covariance matrix NOT Hermitian! Unphysical."        
        w , V = sp_linalg.eigh( C )
        for m in range( w.shape[0] ) :
            w[m] = np.real( w[m] )
            if np.abs(w[m]) / np.max(w) < 1e-10 :
                w[m] = 0
            if w[m] < 0 :
                print 'Negative eigenvalue! Simulating unpysical signal...'

        ntilde_p[ :,k ] =  np.conj( np.sqrt( N / (2*stime) ) * np.dot( V , np.dot( np.sqrt( np.diag( w ) ) , zs[ :,k ] ) ) )
    
    zerofill = np.zeros( ( Nts , 1 ) )
    if N % 2 == 0 :
        ntilde = np.concatenate( ( zerofill , ntilde_p , zerofill , np.conj(np.fliplr(ntilde_p)) ) , axis = 1 )
    else :
        ntilde = np.concatenate( ( zerofill , ntilde_p , np.conj(np.fliplr(ntilde_p)) ) , axis = 1 )
    n = np.real( sp.ifft( ntilde , axis = 1 ) )
    return t , n




def divide_days_in_batches( days , Nb ) :
    """
    Divide the days into batches so they can run in parallel
    INPUT:
    days --- list containing the days. For example: [1,3] for day#1 and day#3
    Nb --- number of batches into which to divide the days
    OUTPUT:
    days_batches --- list containing the days for each batch. For example, in\
    [ [1,3] , [9] ], batch#1 contains day#1 and day#3 and batch#2 contains day#9.
    """
    Ndays = len( days )
    if Ndays % Nb == 0 :
        Ndb = Ndays / Nb
        days_batches = [ days[ b*Ndb : (b+1)*Ndb  ]  for b in range( Nb ) ]
    elif Ndays % Nb > 0 :
        Ndb = Ndays / ( Nb-1 )
        Ndbl = Ndays % ( Nb-1 )
        days_batches = [ days[ b*Ndb : (b+1)*Ndb ] for b in range( Nb-1 ) ] + [ days[ (Nb-1)*Ndb : ] ]
    else :
        raise Exception , "Both the number of days and number of batches have to be postivie integer!"
    return days_batches

def divide_days_in_batches_reverse_pairup( days , Ntogroup ) :
    """
    Divide a list of days in half.  Reverse one half.  Pair the two halves. \n
    Then split the pairs in groups of same number
    INPUT:
    days --- list of days
    Ntogroup --- number of pairs to group together (into a batch)
    OUTPUT
    days_batches --- list whose each element is a list of days for a batch
    """
    if len( days ) % 2 == 0 :
        days_batches = []
    else :
        days_batches = [ [ days.pop() ] ]
    Npairs = len( days ) / 2
    half1st = days[ : Npairs ] ; half2nd = days[ Npairs: ]
    half2nd.reverse()
    pairs = zip( half1st , half2nd )
    pairs = [ list( pair ) for pair in pairs ]
    while True :
        if len( pairs ) >= Ntogroup :
            days_batches += [ list( flatten( [ pairs.pop()
                                               for g in range( Ntogroup ) ] ) ) ]
        elif Ntogroup > len( pairs ) > 0 :
            days_batches += [ list( flatten( pairs ) ) ]
            break
        else :
            break
    return days_batches


def window_and_join( lts , rts , tail=None ) :
    """
    Window two arrays, overlap them by 50% and then add them together.
    INPUT:
    lts --- left array
    rts --- right array
    tail --- array (which maybe output from previous call to the function) to overlap\n
    with the first half of the left array
    OUTPUT:
    ts --- the joined array, made up of the first half of the left array \n
    (to which an input tail may or may not have been added) and the second half\n
    of the left array summed with the first half of the right array 
    tail --- the second half of the right array
    NOTE: All arrays mentioned above are windowed first
    """
    if len( lts.shape ) != 1 or len( rts.shape ) !=1 :
        raise InputError , 'lts and rts must both be 1D arrays!'
    if lts.shape[0] != rts.shape[0] :
        raise InputError , 'lts and rts must have the same length!'
    N = lts.shape[0]
    if N % 2 == 0 :
        midpt = N / 2
    else :
        midpt = ( N - 1 ) / 2
    sinwin = np.sin( np.pi/(N-1) * np.arange( N ) )
    if tail == None :
        ltsw = np.concatenate( ( lts[ :midpt ] , sinwin[ midpt: ]*lts[ midpt: ] ) )
        rtsw = sinwin * rts
        ts = np.concatenate( ( ltsw[ :midpt ] , ltsw[ midpt: ] + rtsw[ :N-midpt ] ) )
        tail = rtsw[ N-midpt: ]
    else :
        if tail.shape != lts[ :midpt ].shape :
            raise InputError, 'The tail must have shape %d!' % lts[ :midpt ].shape[0]
        ltsw = sinwin * lts ; rtsw = sinwin * rts
        ts = np.concatenate( ( tail + ltsw[ :midpt ] ,
                               ltsw[ midpt: ] + rtsw[ :N-midpt ] ) )
        tail = rtsw[ N-midpt: ]
    return ts , tail


def flatten( nested ) :
    try:
        try: nested + ''
        except TypeError: pass
        else: raise TypeError
        for sublist in nested:
            for element in flatten( sublist ) :
                yield element
    except TypeError:
        yield nested


def get_CovarMatrix( Nvar , Ppath , GWSpectralSlope=-3, H0=1.0, *orfpaths ) :
    """
    Returns Nvar x Nvar x Nf covariance matrix from convolving ORF with Plm
    It computes all independent P_{IJ}(f) = H(f)\gamma_{\alpha}^{IJ}(f)P_{\alpha}, and then forms the whole covariance matrix
    INPUT:
    Nvar --- number of variables of covariance matrix (if only A and E, Nvar=2)
    Ppath --- file path to the SkyMap object describing the P_{lm}s
    GWSpectralSlope --- gravitational wave background spectral slope
    H0 --- normalisation constant for H(f)
    *orfpaths --- file paths to the tdiORF_SpHs of the independent cross-correlations
    OUTPUT:
    f --- frequencies ( Nfx1 numpy array )
    comatrix --- covariance matrix ( Nvar x Nvar x  Nf numpy array )
    """
    Nindies = int( Nvar * ( Nvar + 1 ) / 2 )
    if Nindies != len( orfpaths ) :
        raise Exception , 'The covariance matrix %d x %d.  Please make sure that there are %d orfpaths' % ( Nvar , Nvar , Nindies )
    orfs = [ AS.OrfMultipleMoments( orfpath ) for orfpath in orfpaths ]
    f0_df_Nfs = [ ( orf.f.Offset1 , orf.f.Cadence1 , orf.f.data.shape[0] ) for orf in orfs ]
    if not f0_df_Nfs.count( f0_df_Nfs[0] ) == len( f0_df_Nfs ) :
        raise Exception , 'The orfs input do not have the same frequencies.  Please make sure they do!'
    else :
        f = orfs[0].f.data
    file = open( Ppath , 'rb' ) ; skymap = cpkl.load( file ) ; file.close
    QPIJs = [ AS.Convolve( orf , skymap , GWSpectralSlope, H0 ) for orf in orfs ]
    PIJs = [ QPIJ.data for QPIJ in QPIJs ] ; PIJs.reverse()
    comatrix = np.zeros( ( Nvar , Nvar , f.shape[0] ) , dtype = complex )
    for k in range( Nvar ) :
        comatrix[ k , k: , : ] = np.array( [ PIJs.pop() for j in range( Nvar - k ) ] )
        comatrix[ k , :k , : ] = np.copy( np.conj( comatrix[ :k , k , : ] ) )
    return f , comatrix
