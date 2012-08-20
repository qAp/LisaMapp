import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import AnisotropySearch as AS
import myUsefuls as mufls
import scipy  as sp





"""
Test routines for signal/noise simulation.  
"""



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
    """ ###### This is where you define SpHs of real and imaginary parts of ORF. Note that the value returned is for the SpH defined in pyspharm."""  
    indxp = AS.getMLvec( lmax , m='p' )
    SpHimag = np.zeros( ( len(indxp) , f.shape[0] ) , dtype=complex )
    if IJ in [ 'AA' , 'EE' , 'AE' , 'EA' ] :
        SpHreal = np.ones( ( len(indxp) , f.shape[0] ) , dtype=complex ) * 288. / np.sqrt( 2*np.pi )
    elif IJ in [ 'AT','ET','TA','TE','TT' ] :
        SpHreal = np.ones( ( len(indxp) , f.shape[0] ) , dtype=complex ) * 288. / np.sqrt( 2*np.pi )
    """ ###### """
    orfdict = {'OrfMultipleMoments':{ 'ntrunc':lmax , 'f':f , 'Antenna':IJ , 'real': SpHreal , 'imag': SpHimag } }
    orfdir = os.path.dirname( orfpath )
    if orfdir == '' :
        pass
    elif orfdir not in glob.glob( orfdir ) :
        os.system( 'mkdir -p %s' % orfdir )
    file = open( orfpath , 'wb' ) ; cpkl.dump( orfdict , file , -1 ) ; file.close()
    print " Orf multipole moments saved in %s " % orfpath
    return


    
def simulate_AETnoise_from_arbitrary_SpH( duration , stime , t0 ,
                                          GWSpectralSlope , Ppath , lmax ,
                                          seed , N_previous_draws ,
                                          Nvar , compute_ORF_SpHs ,
                                          *orfpaths ) :
    """
    Simulates noise(s) in the freq-domain from SpHs of ORF and SkyMap.
    Partly hardwired to simulate for A,E and T.
    tdiORF_SpHs' IJs are set from AET explicitly, so they can be used for analysis too.
    INPUT:
    duration -- duration of the time-series
    stime -- sampling time 
    t0 --- initial time
    GWSpectralSlope --- GW background spectral slope
    Ppath --- file path to SkyMap of P_{lm}s
    lmax --- maximum degree l in SpH considered
    seed --- seed for the random number generator: 'random' or a positive integer
    N_previous_draws --- numuber of random number draws to discard first
    Nvar --- number of random variables (or time-series)
    compute_ORF_SpHs --- True (to compute tdiORF_SpHs) or False (to load from disk)
    *orfpaths --- file paths to save or to look for tdiORF_SpHs
    OUTPUT:
    t --- time labels of the time-series (N x 1 numpy array) 
    n --- time-series (Nvar x N numpy array)
    """
    Nvar = 3 #3 variables: A, E and T
    IJs = [ 'AA' , 'AE' , 'AT' , 'EE' , 'ET' , 'TT' ] #Explicitly label tdiORF_SpHs (Not really needed here, but this way they can be used for the analysis pipeline.)
    freqdict = mufls.get_freqs_from_duration_and_stime( stime , duration )
    if compute_ORF_SpHs :
        Nindies = int( Nvar*(Nvar+1) / 2 )
        if len( orfpaths ) < Nindies :
            raise InputError , 'You must enter at least %d orfpaths!' % Nindies
        [ make_arbitrary_tdiORF_SpHs( orfpaths[k] , freqdict['f'] , IJ=IJs[k] , lmax=lmax )
          for k in range( Nindies ) ]
        f , comatrix = mufls.get_CovarMatrix( Nvar , Ppath , GWSpectralSlope , *orfpaths )
    else :
        f , comatrix = mufls.get_CovarMatrix( Nvar , Ppath , GWSpectralSlope , *orfpaths )
        if not np.allclose( f , freqdict['f'] ) :
            raise InputError , 'Frequencies of loaded covariance matrix do not match those of the fft.  Try computing the tdiORF_SpHs first with compute_ORF_SpHs=True'
    t , n = mufls.get_noise_freq_domain_CovarMatrix( comatrix , freqdict['df'] , t0 ,
                                                     freqdict['parityN'] ,
                                                     seed=seed ,
                                                     N_previous_draws=N_previous_draws )    
    return t , n


def simulateColorNoise( t0 , deltaT , N , Abar, fR, alpha, seed ) :
    np.random.seed( seed )
    # discrete times
    T = N * deltaT
    t = t0 + deltaT * np.arange( N )
    if N % 2 == 0 :
        numFreqs = N / 2 - 1
    else :
        numFreqs = ( N - 1 ) / 2
    # discrete positive frequencies
    deltaF = 1 / ( N * deltaT )
    f = deltaF + deltaF * np.arange( numFreqs )
    # normalisation factor (discrete freqs > 0)
    norm = np.sqrt(N/(2*deltaT)) * np.sqrt(1/(12*np.pi**2) * Abar**2 * fR**(-2*alpha) * f**(2*alpha-3))

    ys = []
    # need to simulate multiple time-series to avoid periodicity
    for k in range( 3 ) :
        # construct real and imaginary parts, with random phases
        re = (norm/np.sqrt(2)) * np.ones( ( numFreqs, ) )
        im = (norm/np.sqrt(2)) * np.ones( ( numFreqs, ) )
#        re = (norm/np.sqrt(2)) * np.random.standard_normal( ( numFreqs, ) )
#        im = (norm/np.sqrt(2)) * np.random.standard_normal( ( numFreqs, ) )
        z  = re + 1j * im
        # set DC and Nyquist = 0, then add negative freq parts in proper order
        if N % 2 == 0 :
            # note that most negative frequency is -f_Nyquist when N=even
            xtilde = np.array( [0] + list(z) + [0] + list(np.flipud(np.conj(z))) )
        else :
            # no Nyquist frequency when N=odd
            xtilde = np.array( [0] + list(z) + list(np.flipud(np.conj(z))) )
        # fourier transform back to time domain
        temp = sp.ifft( xtilde )
        # take real part (imag part = 0 to round-off)
        ys += [ np.real( temp ) ]

    # splice together the data using sinusoids of twice the period
    w = np.zeros( (N,) )
    for k in range( N ) :
        w[k] = np.sin( ( np.pi * k ) / N )

    y1 = w * ys[0]
    y2 = w * ys[1]
    y3 = w * ys[2]

    z1 = np.array( list( y1[N/2:] ) + list( np.zeros( (N/2,) ) ) )
    z2 = y2
    z3 = np.array( list( np.zeros( (N/2,) ) ) + list( y3[:N/2] ) )
    xdata = z1 + z2 + z3
    return t , xdata

