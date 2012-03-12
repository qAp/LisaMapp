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
    f = df * np.range( 1 , Nf + 1 )
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



def make_arbitrary_tdiORF_SpHs( orfpath , f , IJ='AA' , lmax=0 ) :
    """
    This lets you set multipole moments of overlap-reduction function arbitrarily, and save them in the form\n
    {'OrfMultipleMoments':\n
    { 'ntrunc': lmax , 'f': f , 'Antenna': IJ , 'real': SpHreal , 'imag': SpHimag } }\n
    which is loadable by AnisotropySearch.OrfMultipleMoments()\n

    INPUT:
    orfpath --- file path to save the orfs\n
    f --- frequency/frequencies at which ORF will be computed
    IJ --- IJ of ORF
    lmax --- Maximum degree l for the SpHs
    OUTPUT:
    None
    """
    """ This is where you define SpHs of real and imaginary parts of ORF. Note that the the value returned is for the SpH defined in pyspharm.  As it is, all multipole moments will be the same function of frequency. (The input IJ has no effect on this!)"""
    def arbitrary_SpHreal_ML( f ) :
        g00 = 9
        SpHreal_ML = np.ones( f.shape ) * ( g00/np.sqrt(2*np.pi) )
        return SpHreal_ML
    def arbitrary_SpHimag_ML( f ) :
        SpHimag_ML = np.zeros( f.shape )
        return SpHimag_ML
    indxp = AS.getMLvec( options.lmax , 'p' )
    SpHreal = np.ones( ( len( indxp ) , options.Nf ) ) * arbitrary_SpHreal_ML( f )
    SpHimag = np.ones( ( len( indxp ) , options.Nf ) ) * arbitrary_SpHimag_ML( f )
    orfdict = {'OrfMultipleMoments': { 'ntrunc': options.lmax , 'f': f , 'Antenna': options.IJ , 'real': SpHreal , 'imag': SpHimag } }
    orfdir = os.path.dirname( orfpath )
    if orfdir not in glob.glob( orfdir ) :
        os.system( 'mkdir -p %s' % orfdir )
    file = open( orfpath , 'wb' ) ; cpkl.dump( orfdict , file , -1 ) ; file.close()
    print " Orf multipole moments saved in %s " % orfpath
    return


def get_CovarMatrix( Nvar , Ppath , GWSpectralSlope=-3 , *orfpaths ) :
    """
    Returns Nvar x Nvar x Nf covariance matrix from convolving ORF with Plm
    It computes all independent P_{IJ}(f) = H(f)\gamma_{\alpha}^{IJ}(f)P_{\alpha}, and then forms the whole covariance matrix
    INPUT:
    Nvar --- number of variables of covariance matrix (if only A and E, Nvar=2)
    Ppath --- file path to the SkyMap object describing the P_{lm}s
    GWSpectralSlope --- gravitational wave background spectral slope
    *orfpaths --- file paths to the tdiORF_SpHs of the independent cross-correlations
    OUTPUT:
    f --- frequencies ( Nfx1 numpy array )
    comatrix --- covariance matrix ( Nvar x Nvar x  Nf numpy array )
    """
    Nindies = int( Nvar * ( Nvar + 1 ) / 2 )
    if Nindies != len( orfpaths ) :
        raise Exception , 'The covariance matrix %d x %d.  Please make sure that there are %d orfpaths' % ( Nvar , Nvar , Nindies )
    orfs = [ OrfMultipleMoments( orfpath ) for orfpath in orfpaths ]
    f0_df_Nfs = [ ( orf.f.Offset1 , orf.f.Cadence1 , orf.f.data.shape[0] ) for orf in orfs ]
    if not f0_df_Nfs.count( f0_df_Nfs[0] ) == len( f0_df_Nfs ) :
        raise Exception , 'The orfs input do not have the same frequencies.  Please make sure they do!'
    else :
        f = orfs[0].f.data
    file = open( Ppath , 'rb' ) ; skymap = cpkl.load( file ) ; file.close
    QPIJs = [ AS.Convolve( orf , skymap , options.GWslope ) for orf in orfs ]
    PIJs = [ QPIJ.data for QPIJ in QPIJs ] ; PIJs.reverse()
    comatrix = np.zeros( ( Nvar , Nvar , f.shape[0] ) , dtype = complex )
    for k in range( Nvar ) :
        comatrix[ k , k: , : ] = np.array( [ PIJs.pop() for j in range( Nvar - k ) ] )
        comatrix[ k , :k , : ] = np.copy( np.conj( comatrix[ :k , k , : ] ) )
    return f , comatrix

    
def simulate_noise_from_arbitrary_SpH( duration , stime , t0 ,
                                       GWSpectralSlope , Ppath ,
                                       seed='random' , N_previous_draws = 0 ,
                                       Nvar=0 , compute_ORF_SpHs = False ,
                                       *orfpaths ) :
    freqdict = get_freqs_from_duration_and_stime( stime , duration )
    if compute_ORF_SpHs :
        pass
    else :

        
    return 
