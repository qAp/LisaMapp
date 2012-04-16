import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import AnisotropySearch as AS
import myUsefuls as mufls






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
    """ This is where you define SpHs of real and imaginary parts of ORF. Note that the the value returned is for the SpH defined in pyspharm.  As it is, all multipole moments will be the same function of frequency. (The input IJ has no effect on this!)"""
    def arbitrary_SpHreal_ML( f ) :
        g00 = 144.
        SpHreal_ML = np.ones( f.shape ) * ( g00/np.sqrt(2*np.pi) )
        return SpHreal_ML
    def arbitrary_SpHimag_ML( f ) :
        SpHimag_ML = np.zeros( f.shape )
        return SpHimag_ML
    indxp = AS.getMLvec( lmax , 'p' )
    SpHreal = np.ones( ( len( indxp ) , f.shape[0] ) ) * arbitrary_SpHreal_ML( f )
    SpHimag = np.ones( ( len( indxp ) , f.shape[0] ) ) * arbitrary_SpHimag_ML( f )
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
