
#import synthlisa
import lisaxml
import numpy
import LISAresponse
import pylab
import sys
import glob
import os
import re
import cPickle
import matplotlib.pyplot as pyplot
import Utilities4 as U4


def sumTerms(x,jlow,jhigh):
    """Return the sums of elements in x between indices jlow and jhigh.

    Positional parameters:
    x -- numpy array 
    jlow -- numpy array containing lower indices
    jhigh -- numpy array containing upper indices

    Output:
    y -- numpy array of same length as jlow and jhigh; each element in the ouput array corresponds to an index in jlow and in jhigh, and is the sum of elements in x between these two indices.

    """
    Ny = len(jlow)
    y = numpy.zeros(Ny,complex)
    for k in range(Ny):
        y[k] = sum(x[jlow[k]:jhigh[k]])
    return y

def coarsegrain(data,flowx,deltafx,flowy,deltafy,Ny):
    """Return coarsegrained data

    Positional parameters:
    data    -- numpy array (either 1D or 2D)
    flowx   -- data's offset
    deltafx -- data's cadence
    flowy   -- desired offset
    deltafy -- desired cadence
    Ny      -- desired length

    Output:
    input array coarsgrained to the desired offset, cadence and length.  

    Note:
    The function always coarsegrains along the 2nd dimension of the array.  If data is 1D, it will first be converted into a 2D array with one row. However, it will be returned as a coarsegrained 1D array.  

    """
    if len( data.shape ) == 1 :
        data = numpy.array( [ data ] )
        datawas1D = True
    elif len( data.shape ) == 2 :
        datawas1D = False
        pass
    else:
        print "The array to be coarse-grained can only be 1-D or 2-D. "
        raise

    Nx = data.shape[1]
    if  (deltafx <= 0) or (deltafy <= 0) or (Ny <= 0):
        print 'bad input argument'
        raise ValueError 
    if  deltafy < deltafx :
        print 'deltaf coarse grain < deltaf fine grain'
        raise ValueError 
    if "%.10f" % (flowy - 0.5*deltafy) < "%.10f" % (flowx - 0.5*deltafx):
        print 'desired coarse grained start frequency is too low'
        raise ValueError 
    fhighx = flowx + (Nx-1)*deltafx
    fhighy = flowy + (Ny-1)*deltafy
    if "%.10f" % (fhighy + 0.5*deltafy) > "%.10f" % (fhighx + 0.5*deltafx):
        print 'desired coarse grained stop frequency is too high'
        raise ValueError 
    i = numpy.arange(Ny)
    jlow = numpy.intp(
        numpy.ceil(
        numpy.round( ( flowy + (i-0.5)*deltafy  -  flowx + 0.5*deltafx ) / deltafx , 10 )
        ) - 1
        )
    jhigh = numpy.intp(
        numpy.floor(
        numpy.round( ( flowy + (i+0.5)*deltafy - flowx + 0.5*deltafx  ) / deltafx , 10 )
        )
        )
    index1 = jlow[0]
    index2 = jhigh[-1]
    fraclow  = ( flowx+(jlow+0.5)*deltafx - flowy-(i-0.5)*deltafy  ) / deltafx
    frachigh = ( flowy+(i+0.5)*deltafy - flowx-(jhigh-0.5)*deltafx ) / deltafx
    frac1 = fraclow[0]
    frac2 = frachigh[-1]
    jtemp = jlow + 1

    coarsedata = numpy.zeros( (data.shape[0],Ny) , complex )
    for lm in range( data.shape[0] ):
        midsum = sumTerms( data[lm,:] , jtemp , jhigh )
        ya = (float(deltafx)/deltafy) * ( data[lm,:][jlow[:-1]] *fraclow[:-1]  + 
                                          data[lm,:][jhigh[:-1]]*frachigh[:-1] + 
                                          midsum[:-1]
                                          )
        if ( jhigh[-1] > Nx-1 ):
            yb = (float(deltafx)/deltafy) * ( data[lm,:][jlow[-1]]*fraclow[-1] + 
                                              midsum[-1]
                                              )
        else:
            yb = (float(deltafx)/deltafy) * ( data[lm,:][jlow[-1]] *fraclow[-1] +
                                              data[lm,:][jhigh[-1]]*frachigh[-1] + 
                                              midsum[-1]
                                              )

        coarsedata[lm,:] = numpy.array( list(ya) + [yb] )
    if datawas1D :
        coarsedata = coarsedata[0]
    return coarsedata


class Coarsable(object):
    """Coarsable is for making a data array of 1D or 2D coarsegrainable.

    Attribute:
    data -- 1D/2D numpy array containing data
    Offset0  -- offset along data's first dimension
    Cadence0 -- cadence along data's first dimension
    Offset1  -- offset along data's second dimension
    Cadence1 -- cadence along data's second dimension
    frame -- dictionary containing the offsets and cadences

    Methods:
    coarsegrain -- convert data to correpond to coarser cadences
    
    """

    def __init__(self,
                 data,
                 Offset0=None,Cadence0=None,Offset1=None,Cadence1=None):
        """Create an instance of the class Coarsable

        Positional parameters:
        data -- 1D/2D numpy array 

        Keyword parameters:
        Offset0  -- offset along data's first dimension (default None)
        Cadence0 -- cadence along data's first dimension (default None)
        Offset1  -- offset along data's second dimension (default None)
        Cadence1 -- cadence along data's second dimension (default None)

        Output:
        a coarsable object

        """

        self.data = data
        self.Offset0  = Offset0
        self.Cadence0 = Cadence0
        self.Offset1  = Offset1
        self.Cadence1 = Cadence1
        self.frame = { 'Offset0':self.Offset0 , 'Cadence0':self.Cadence0 ,
                       'Offset1':self.Offset1 , 'Cadence1':self.Cadence1 }
        return
    def coarsegrain(self,
                    Offset0=None,Cadence0=None,N0=None,
                    Offset1=None,Cadence1=None,N1=None):
        """Return coarsegrained data 

        Keyword parameters:
        Offset0  -- desired offset of data's first dimension (default None)
        Cadence0 -- desired cadence of data's first dimension (default None)
        N0       -- desired length of data's first dimension (default None)
        Offset1  -- desired offset of data's second dimension (default None)
        Cadence1 -- desired cadence of data's second dimension (default None)
        N1       -- desired length of data's second dimension (default None)

        Output:
        coarsable object with new offsets and cadences, whose data is coarsegrained accordingly.

        Note:
        Coarsegraining is possible along both sides of the data array. Their order ought not to matter.

        """

        if None in [Offset0,Cadence0,N0] and None in [Offset1,Cadence1,N1]:
            print "Nothing to do. Try specifying coarse bins' starting offset, cadence and length."
            coarsedata = Coarsable( self.data , **self.frame )
            return coarsedata
        elif None in [Offset0,Cadence0,N0] and None not in [Offset1,Cadence1,N1]:
            if None in [self.Offset1,self.Cadence1]:
                print "Cannot coarsegrain side 1.  No offset and cadence assigned to it."
                coarsedata = Coarsable( self.data , **self.frame )
            else:
                try: 
                    print "Trying to coarsegrain side 1..."
                    coarsedata = coarsegrain(self.data,self.Offset1,self.Cadence1,Offset1,Cadence1,N1)
                    print "Side 1 coarsegrained."
                    coarsedata = Coarsable( coarsedata , 
                                            Offset0=self.Offset0 , Cadence0=self.Cadence0 , 
                                            Offset1=Offset1 , Cadence1=Cadence1 )
                except ValueError:
                    coarsedata = Coarsable( self.data , **self.frame )
                    print "Side 1 not coarsegrained."
            return coarsedata
        elif None in [Offset1,Cadence1,N1] and None not in [Offset0,Cadence0,N0]:
            if None in [self.Offset0,self.Cadence0]:
                print "Cannot coarsegrain side 0.  No offset and cadence assigned to it."
                coarsedata = Coarsable( self.data , **self.frame )
            else:
                try:
                    print "Trying to coarsegrain side 0..."
                    data       = numpy.transpose(self.data)
                    coarsedata = coarsegrain(     data,self.Offset0,self.Cadence0,Offset0,Cadence0,N0)
                    print "Side 0 coarsegrained."
                    coarsedata = Coarsable( numpy.transpose(coarsedata) , 
                                            Offset0=Offset0 , Cadence0=Cadence0 ,
                                            Offset1=self.Offset1 , Cadence1=self.Cadence1 )
                except ValueError:
                    coarsedata = Coarsable( self.data , **self.frame )
                    print "Side 0 not coarsegrained."
            return coarsedata
        else:
            coarse0data = self.coarsegrain( Offset0=Offset0 , Cadence0=Cadence0 , N0=N0 ,
                                            Offset1=None , Cadence1=None , N1=None )
            coarsedata = coarse0data.coarsegrain( Offset0=None , Cadence0=None , N0=None ,
                                                  Offset1=Offset1 , Cadence1=Cadence1 , N1=N1 ) 
            return coarsedata
#        else:
#            if None in [self.Offset0,self.Cadence0]:
#                print "Cannot coarsegrain side 0. No offset and cadence assigned to it."
#                coarse0data = Coarsable( self.data , **self.frame )
#            else:
#                try:
#                    print "Trying to coarsegrain side 0..."
#                    data        = numpy.transpose(self.data)
#                    coarse0data = coarsegrain(data,self.Offset0,self.Cadence0,Offset0,Cadence0,N0)
#                    print "Side 0 coarsegrained."
#                    coarse0data = Coarsable( numpy.transpose(coarse0data) , Offset0=Offset0 , Cadence0=Cadence0 )
#                except ValueError:
#                    coarse0data = Coarsable( self.data , **self.frame )
#                    print "Side 0 not coarsegrained."
#            if None in [self.Offset1,self.Cadence1]:
#                print "Cannot coarsegrain side 1. No offset and cadence assigned to it."
#                coarsedata = coarse0data
#            else:
#                try:
#                    print "Trying to coarsegrain side 1..."
#                    coarsedata = coarsegrain(coarse0data.data,self.Offset1,self.Cadence1,Offset1,Cadence1,N1)
#                    print "Side 1 coarsegrained."
#                    coarsedata = Coarsable( coarsedata ,
#                                            Offset0=coarse0data.Offset0 , Cadence0=coarse0data.Cadence0 ,
#                                            Offset1=Offset1 , Cadence1=Cadence1 )
#                except ValueError:
#                    coarsedata = coarse0data
#                    print "Side 1 not coarsegrained."
#        return coarsedata



def coarsefrequency(*fcoarsables):
    """Return the coarsest frequencies in the common frequency range.
    
    Positional parameter:
    *fcoarsables -- coarsable object/objects of a 1D numpy array of frequencies.

    Output:
    Coarsable object of a 1D numpy array of frequencies. 
    The cadence is the largest among the input set of frequencies.
    The range is the widest and common range shared by the input set of frequencies.
    
    """

    offsets = [ f.Cadence1 for f in fcoarsables ] 
    offsets_coarsables = zip( offsets , fcoarsables )
    offsets_coarsables.sort()
    offsets_coarsables.reverse()
    coarsables = [ f[1] for f in offsets_coarsables ]
    F , fineones = coarsables[0] , coarsables[1:]

    leds  = [ f.data[0]  - f.Cadence1/2. for f in fineones ]
    leds_fineones = zip( leds , fineones )
    leds_fineones.sort()
    highest_led , lowone  = leds_fineones[-1][0] , leds_fineones[-1][-1]

    heds = [ f.data[-1] + f.Cadence1/2. for f in fineones ]
    heds_fineones = zip( heds , fineones )
    heds_fineones.sort()
    lowest_hed , highone = heds_fineones[0][0] , heds_fineones[0][-1]
        
    if F.data[0]-F.Cadence1/2. >= highest_led :
        Nlow = 0
    else:
        Nlow = numpy.ceil( (highest_led - (F.data[0]-F.Cadence1/2.)) / F.Cadence1 ) 

    if F.data[-1]+F.Cadence1/2. <= lowest_hed :
        data = F.data[ Nlow : ]
    else:
        Nhigh = numpy.floor( (lowest_hed - (F.data[0]+F.Cadence1/2.)) / F.Cadence1 ) 
        data = F.data[ Nlow : Nhigh + 1 ]

    fcoarse = Coarsable( data , 
                         Offset1 = data[0] , Cadence1 = F.Cadence1 )
    return fcoarse


def tsXML_to_tsDICT( tspath ):
    inxml = lisaxml.readXML( tspath )
    tdiobs = inxml.TDIData[0]

    ts_scale = { 'Offset1'  : tdiobs.TimeSeries.TimeOffset ,
                 'Cadence1' : tdiobs.TimeSeries.Cadence
                }
    
    tsdict = {}
    tsdict['t'] = Coarsable( tdiobs.t , **ts_scale )
    tsdict['A'] = Coarsable( tdiobs.A , **ts_scale )
    tsdict['E'] = Coarsable( tdiobs.E , **ts_scale )
    tsdict['T'] = Coarsable( tdiobs.T , **ts_scale )
    return tsdict


class TimeSeries(object):
    """TimeSeries class handles time-series in the analysis

    Attributes:
    t -- time 
    A -- time-series A
    E -- time-series E
    T -- time-series T (all the above 4 attributes are Coarsables)

    Methods:
    psd -- Return the PSD(power spectral density) of the time-series
    csd -- Return the CSD(cross spectral density) of the time-series
    winfft -- Window the time-series, and then return the fft
    CrossSpectraData -- Return cross spectral data

    """

    def __init__(self,tsdict):
        """Create an instance of the class TimeSeries.

        Positional parameter:
        tsdict -- dictionary containing t,A,E,T

        OUTPUT:
        A TimeSeries object containing time-series info

        """
        self.t = tsdict['t']
        self.A = tsdict['A']
        self.E = tsdict['E']
        self.T = tsdict['T']
        return
    def psd(self,xstr,nfft=256,noverlap=0,coarsable=False):
        """Return the PSD(power spectral density) of the time-series

        Positional parameter:
        xstr -- string naming the observable, either 'A', 'E', or 'T'

        Keyword parameters:
        nfft -- length of the fast fourier transform (fft) (default 256)
        noverlap -- number of data points overlapping between segments (default 0)
        coarsable -- return the PSD as a Coarsable object (default False)
        
        OUTPUT:
        f -- frequencies 
        P -- PSD 

        Note:
        The PSD is calculated using pylab.psd() with an additional factor of dt in the normalisation.
        This gives the correct dimension.

        """
        x = getattr(self,xstr)
        fs = 1. / x.Cadence1
        P , f = pylab.psd( x.data , nfft , fs , noverlap=noverlap , scale_by_freq=False )
        P = P / fs
        if coarsable:
            scale = {'Offset1' : f[0] ,
                     'Cadence1': f[1] - f[0] }
            f = Coarsable( f , **scale )
            P = Coarsable( P , **scale )
        return f , P
    def csd(self,xstr,ystr,nfft=256,noverlap=0,coarsable=False):
        """Return CSD(cross spectral density) of the time-series.

        Positional parameters:
        xstr -- string naming one observable, either 'A', 'E' or 'T'
        ystr -- string naming one observable, either 'A', 'E' or 'T'

        Keyword parameters:
        nfft -- length of fast fourier transform(fft) (default 256)
        noverlap -- number of data points overlapping between segments (default 0)
        coarsable -- return the CSD as a Coarsable object (default False)

        OUTPUT:
        f -- frequencies
        P -- CSD

        Note:
        The CPSD is calculated using pylab.csd() with an additional factor of dt in the normalisation.
        If xstr='A' and ystr='E', then the CSD between A and E is returned.

        """
        x = getattr(self,xstr)
        y = getattr(self,ystr)
        fs = 1. / x.Cadence1
        P , f = pylab.csd( x.data , y.data , nfft , fs , noverlap=noverlap )
        P = P / fs
        if coarsable:
            scale = {'Offset1' : f[0] ,
                     'Cadence1': f[1] - f[0] }
            f = Coarsable( f , **scale )
            P = Coarsable( P , **scale )
        return f , P
    def winfft(self,xstr,window,fftlength):
        """Window the time-series, and then return the fft

        Positional parameters:
        xstr -- string naming the observable, either 'A', 'E' or 'T'
        window -- numpy array containing the amplitude of the window
        fftlength -- length of fast fourier transform(fft)

        OUTPUT:
        f -- frequencies in a numpy array
        winfftx -- Windowed and then fft'ed data in a numpy array

        Note:
        The length of fft has to be longer than the data's length, where the data is zero padded before it is fft'ed.
        The window's length must match the data's length.

        """
        x = getattr( self , xstr )
        if ( len(x.data) != len(window) ):
            raise Exception, 'size mismatch, length of data != that of window'
        if ( fftlength < len(x.data) ):
            raise ValueError, 'fftlength < length of data' 
        z = numpy.zeros( fftlength - len(x.data) )
        xbar = pylab.fft( numpy.array( list(x.data*window) + list(z) ) )
        winfftx = xbar[ 0:numpy.floor( fftlength/2+1 ) ]*x.Cadence1
        df = 1. / (x.Cadence1*fftlength)
        f =  df*numpy.arange( len(winfftx) )
        return f , winfftx
    def CrossSpectraData(self,xstr,ystr,window,fftlength,coarsable=False):
        """Return cross spectral data 

        Positional parameters:
        xstr -- observable, either 'A', 'E' or 'T'
        ystr -- observable, either 'A', 'E' or 'T'
        window -- amplitudes of window for the data, in a numpy array
        fftlength -- length of fft

        Keyword parameter:
        coarsable -- return the cross spectral data as a Coarsable object (default False)

        OUTPUT:
        f -- frequencies
        Cxy -- cross spectral data

        """
        x = getattr( self , xstr )
        dataduration = x.Cadence1 * len(x.data)
        f , winfftx = self.winfft( xstr , window , fftlength )
        f , winffty = self.winfft( ystr , window , fftlength )
        Cxy = numpy.conj(winfftx)*winffty * 2 / dataduration
        if coarsable:
            scale = {'Offset1' : f[0] ,
                     'Cadence1': f[1] - f[0] }
            f   = Coarsable( f , **scale )
            Cxy = Coarsable( Cxy , **scale )
        return f , Cxy


class OrfMultipleMoments(object):
    """OrfMultipleMoments handles the overlap-reduction function(orf)'s multiple moments

    Attributes:
    LISAData -- type and initial configuration of LISA
    Sky -- sky object created from spharm.Spharmt, like the grid used to cover the sky
    Antenna -- the pair of detectors used in the cross-correlation
    ntrunc -- lmax
    fscale -- dictionary containing the frequency offset and cadence
    f -- Coarsable of the frequencies
    real -- Coarsable containing the multiple moments of orf's real part
    imag -- Coarsable containing the multiple moments of orf's imag part

    Methods:
    getMultipleMoments -- return multiple moments of orf

    """
    def __init__(self,orfpath):
        """Create an instance of the class OrfMultipleMoments

        Positional parameters:
        orfpath -- absolute/relative path of pickle(.pkl) file containing orf's multiple moments

        OUTPUT:
        An object of class OrfMultipleMoments
        """
        file = open( orfpath , 'rb' )
        orf  = cPickle.load( file )
        file.close()
        self.LISAData = orf['LISAData']
        self.Sky      = orf['Sky']
        self.Antenna  = orf['OrfMultipleMoments']['Antenna']
        self.ntrunc   = orf['OrfMultipleMoments']['ntrunc']
        f    = orf['OrfMultipleMoments']['f']
        real = orf['OrfMultipleMoments']['real']
        imag = orf['OrfMultipleMoments']['imag']
        self.fscale = { 'Offset1':f[0] , 'Cadence1':f[1]-f[0]}
        self.f    = Coarsable( f    , **self.fscale )
        self.real = Coarsable( real , **self.fscale )
        self.imag = Coarsable( imag , **self.fscale )
        return
    def getMultipleMoments(self,msign='pn',lmax=20):
        """Return multiple moments of orf

        Keyword parameters:
        msign -- 'p'  -- return only positive m's multiple moments
              -- 'pn' -- return all multiple moments

        Output:
        g -- Coarsable object containing orf's multiple moments

        """
        sqrt2pi = numpy.sqrt( 2*numpy.pi )
        if msign == 'p' :
            print "Warning: input lmax ineffective for msign='p'"
            p , q = numpy.conj( self.real.data ) , numpy.conj( self.imag.data )
            indxp = LISAresponse.getMLvec( self.ntrunc , 'p' )
            data = numpy.zeros( p.shape , dtype='complex' )
            for i,ml in enumerate( indxp ):
                data[i] = (-1)**ml[0] * sqrt2pi * ( p[i] + 1j*q[i] )
            g = Coarsable( data , **self.fscale )
        elif msign == 'pn' :
            p , q = numpy.conj( self.real.data ) , numpy.conj( self.imag.data )
            indxp  = LISAresponse.getMLvec( self.ntrunc , 'p' )
            indxpn = LISAresponse.getMLvec( self.ntrunc , 'pn' )
            data = numpy.zeros( ( len(indxpn) , p.shape[1] ) , complex )
            for i,ml in enumerate(indxpn):
                m , l = ml[0] , ml[1]
                if m >= 0 :
                    k = indxp.index( ml )
                    data[i,:] = (-1)**m * sqrt2pi * ( p[k,:] + 1j*q[k,:] )
                else:
                    ml = (-m,l)
                    k = indxp.index( ml )
                    data[i,:] =  sqrt2pi * ( numpy.conj( p[k,:] ) + 1j*numpy.conj( q[k,:] ) )
            data = get_lmax_subset_from( data , lmax , axis=0 )
            g = Coarsable( data , **self.fscale )
        else:
            print "Keyword parameter 'msign' must either be 'p' for positive m's or 'pn' for all m's."
            raise
        return g


class XXX(object):
    """XXX is for dirty map estimation

    Attributes:
    ts -- time-series in a TimeSeries object
    psd -- PSD in a dictionary
    orf -- overlap-reduction function in a OrfMultipleMoments object
    alpha -- GW spectral slope
    pair -- detectors in the cross-correlation
    fcoarse -- frequencies used foarsegraining.
    
    Methods:
    getCrossSpectraData -- return cross spectral data
    getCoarseFreqs -- return frequencies for coarsegraining
    getintegrand -- return the integrand in the expression for the dirty map
    getsummand -- return the summand in the expression for the dirty map

    Note:
    The class is designed to evaluate expression () in "[paper/note name]".  The various terms in this expression come from the different input objects in the constructor.  Their frequencies are turned into a common set, and the expression is evaluated at these frequencies.

    """
    def __init__(self, OrfMultipleMoments , PSDDIC , TimeSeries , 
                 GwSpectralSlope , sourcename , day , cIJdir='./' ,
                 XXXdir='./' ):
        """Create an instance of the class XXX
        
        Positional parameters:
        OrfMultipleMoments -- OrfMultipleMoments object
        PSDDIC -- PSD dictionary(these are saved in a .pkl file)
        TimeSeries -- TimeSeries object
        GWSpectralSlope -- GW spectral slope
        
        """
        self.ts      = TimeSeries
        self.psd     = PSDDIC
        self.orf     = OrfMultipleMoments
        self.alpha   = GwSpectralSlope
        self.pair    = self.orf.Antenna
######################################################################
        self.sourcename = sourcename
        self.day = day
        self.cIJdir = cIJdir
        self.XXXdir = XXXdir
######################################################################
        return
    def getCrossSpectraData(self,window=None,fftlength=None):
        """Return cross spectral data

        Keyword parameters:
        window -- amplitude of window in numpy array
        fftlength -- length of fast fourier transform (fft)

        Output:
        f -- frequency
        C -- cross spectral data

        """
        if window == None :
            window = pylab.hanning( len( self.ts.t.data ) )
        if fftlength == None :
            fftlength = len( self.ts.t.data )
        inputs = ( self.pair[0] , self.pair[1] , window , fftlength )
        f , C = self.ts.CrossSpectraData( coarsable=True , *inputs )
        self.fcsdata , self.csdata = f , C
        return f , C
    def getCoarseFreqs(self):
        """Return frequencies for coarsegraining

        Output:
        fcoarse -- frequencies in a Coarsable
        
        Note:
        This takes not input parameters. It also assigns the coarsegrain frequencies as a new attirbute of the instance.
        
        """
        if not hasattr(self,'fcsdata') or not hasattr(self,'csdata'):
            self.getCrossSpectraData()
        fcoarse = coarsefrequency( self.orf.f , 
                                   self.psd['PowerSpectra']['f'] , 
                                   self.fcsdata )
        self.fcoarse = fcoarse
        return fcoarse

    def getintegrand(self,flow=1.4e-4,fhigh=1e-1 , lmax=20 ):
        """Return the integrand in the expression for the dirty map

        Keyword parameter:
        flow -- lower frequency limit (default 1.4e-4Hz)
        fhigh -- higher frequency limit (default 1e-1Hz)
        
        Output:
        f -- frequencies
        integrand -- the integrand

        Note:
        See expression () in "[paper/notes name]"
        
        """
        print 'Calculating XXX...'
        if not hasattr( self , 'fcoarse' ) :
            self.getCoarseFreqs()

        if  fhigh < flow :
            print "fhigh must be >= flow"
            raise 
        elif flow < self.fcoarse.data[0] :
            print "flow must be >= %f" % self.fcoarse.data[0]
            raise 
        elif fhigh > self.fcoarse.data[-1] :
            print "fhigh must be <= %f" % self.fcoarse.data[-1]
            raise 

        glm = self.orf.getMultipleMoments( lmax=lmax )
        pII = self.psd['PowerSpectra'][ 'P'+2*self.pair[0] ]
        pJJ = self.psd['PowerSpectra'][ 'P'+2*self.pair[1] ]
        csdata = self.csdata
        
###################### save csdata , pII , pJJ ########################(2)
#        cIJdict = {}
#        cIJdict['csdata']  = csdata
#        cIJdict['pII']     = pII
#        cIJdict['pJJ']     = pJJ
#        
#        file = open( self.cIJdir + self.sourcename + '_obs-s%03d-cIJ.pkl' % 
#                     self.day , 'wb' )
#        cPickle.dump( cIJdict , file , -1 )
#        file.close()
#########################################################################
        
        inputs = { 'Offset1'  : self.fcoarse.Offset1 ,
                   'Cadence1' : self.fcoarse.Cadence1 ,
                   'N1'       : len( self.fcoarse.data ) }
        glm = glm.coarsegrain( **inputs )
        pII , pJJ  = pII.coarsegrain( **inputs ) , pJJ.coarsegrain( **inputs )
        csdata = csdata.coarsegrain( **inputs )

######################## save csdata , pII , PJJ ####################(1)
        cIJdict = {}
        cIJdict['fcoarse'] = self.fcoarse
        cIJdict['csdata']  = csdata
        cIJdict['pII']     = pII
        cIJdict['pJJ']     = pJJ
        file = open( self.cIJdir + self.sourcename + '_obs-s%03d-cIJ.pkl' %
                     self.day , 'wb' )
        cPickle.dump( cIJdict , file , -1 )
        file.close()
###########################################################################

        H = self.fcoarse.data**self.alpha

########## Estimator containing only +ve frequencies ######################
#        HoPP = H / ( numpy.real( pII.data ) * numpy.real( pJJ.data ) )
#        data = numpy.conj( glm.data ) * HoPP * csdata.data
###########################################################################
########## Estimator containing +ve and -ve frequencies ###################
        datap = numpy.conj( glm.data )
        HoPP = H / ( numpy.real( pII.data ) * numpy.real( pJJ.data ) )

        datan = numpy.zeros( datap.shape , dtype='complex' )
        indxpn = LISAresponse.getMLvec( lmax , 'pn' )
        for i,ml in enumerate( indxpn ) :
            m , l = ml
            k = indxpn.index( ( -m , l ) )
            datan[ i ] = (-1)**m * glm.data[ k ]
                        
        data = ( datap * csdata.data + datan * numpy.conj( csdata.data ) ) * HoPP
###########################################################################
        
        ilow  = numpy.round( (flow - self.fcoarse.Offset1) / 
                             self.fcoarse.Cadence1 )
        ihigh = numpy.round( (fhigh - self.fcoarse.Offset1) / 
                             self.fcoarse.Cadence1 )
        data = data[ : , ilow : ihigh + 1 ]
        fdata = self.fcoarse.data[ ilow : ihigh + 1 ]
        fscale = {'Offset1' : fdata[0] , 
                  'Cadence1': fdata[1] - fdata[0] }
        f         = Coarsable( fdata , **fscale ) 
        integrand = Coarsable( data  , **fscale )
        return f , integrand
    def getsummand(self,flow=1.4e-4,fhigh=1e-1 , lmax=20 ):
        """Return the summand in the expression for the dirty map

        Keyword parameters:
        flow -- lower frequency limit (default 1.4e-4Hz)
        fhigh -- higher frequency limit (default 1e-1Hz)

        Output:
        summand -- the summand

        Note:
        See expression () in "[paper/note name]".  The summand returned here is the daily dirty map.
        
        """
        print 'Calculating XX...'
        f , integrand = self.getintegrand( flow , fhigh , lmax=lmax )
################ Extract XXX (1) ###############################
#        XXXdict = {}
#        XXXdict['f']   = f
#        XXXdict['XXX'] = integrand
#        file = open( self.XXXdir + self.sourcename +
#                     '_obs-s%03d-XXX.pkl' % self.day , 'wb' )
#        cPickle.dump( XXXdict , file , -1 )
#        file.close()
################################################################
        data = numpy.sum( integrand.data , 1 )
        summand = Coarsable( data )
        print 'done'
        return summand


class GGG(object):
    """GGG is for getting the ingredients for the Fisher matrix

    Attributes:
    

    Methods:

    
    """
    def __init__(self, OrfMultipleMoments , PSDDIC , GwSpectralSlope , day ):
        """Create an instance of the class GGG
        
        Positional parameters:
        OrfMultipleMoments -- OrfMultipleMoments object
        PSDDIC -- PSD dictionary (these are saved in .pkl files)
        GwSpectralSlope -- GW spectral slope
        
        """
        self.orf   = OrfMultipleMoments
        self.psd   = PSDDIC
        self.alpha = GwSpectralSlope
        self.pair  = self.orf.Antenna

        self.day = day
        return
    def getCoarseFreqs(self):
        """Return frequencies for coarsegraining 

        Output:
        fcoarse -- frequencies in a Coarsable           

        Note:
        This takes not input parameters. It also assigns the coarsegrain frequencies as a new attirbute of the instance.         

        """
        fcoarse = coarsefrequency( self.orf.f ,
                                   self.psd['PowerSpectra']['f'] )
        self.fcoarse = fcoarse
        return fcoarse

    def getmlmlintegrand( self , ml1 , ml2 , flow=1.4e-4 , fhigh=1e-1 , lmax=20 ):
        """Return the integrand in the expression for the Fisher matrix

        Positional parameters:
        ml1 -- (m,l) tuple
        ml2 -- (m,l) tuple

        Keyword parameters:
        flow -- lower frequency limit
        fhigh -- higher frequency limit

        Output:
        f -- frequencies
        mlmlintegrand -- the integrand

        """
        print 'Calculating GGG...'
        if not hasattr( self , 'fcoarse' ):
            self.getCoarseFreqs()
        if fhigh < flow :
            print "fhigh must be >= flow"
            raise
        elif flow < self.fcoarse.data[0] :
            print "flow must be >= %f" % self.fcoarse.data[0]
            raise
        elif fhigh > self.fcoarse.data[-1] :
            print "fhigh must be <= %f" % self.fcoarse.data[-1]
            raise

        glm = self.orf.getMultipleMoments( lmax=lmax )
        pII = self.psd['PowerSpectra'][ 'P'+2*self.pair[0] ]
        pJJ = self.psd['PowerSpectra'][ 'P'+2*self.pair[1] ]
        inputs = { 'Offset1'  : self.fcoarse.Offset1 ,
                   'Cadence1' : self.fcoarse.Cadence1 ,
                   'N1'       : len( self.fcoarse.data ) }
        glm = glm.coarsegrain( **inputs )
        pII , pJJ = pII.coarsegrain( **inputs ) , pJJ.coarsegrain( **inputs )
        H = self.fcoarse.data**self.alpha

        ml1m , ml2m = ( -ml1[0] , ml1[-1] ) , ( -ml2[0] , ml2[-1] )
        indxpn = LISAresponse.getMLvec( self.orf.ntrunc , 'pn' )
        k1 , k2 = indxpn.index( ml1 ) , indxpn.index( ml2 )
        k3 , k4 = indxpn.index( ml1m ) , indxpn.index( ml2m )
        
        data = ( numpy.conj( glm.data[k1,:] ) * glm.data[k2,:] + (-1)**( ml1[0] + ml2[0] ) * numpy.conj( glm.data[k4,:] ) * glm.data[k3,:] ) * H**2 / numpy.real( pII.data ) / numpy.real( pJJ.data )

        ilow  = int( numpy.round( (flow - self.fcoarse.Offset1) /
                                  self.fcoarse.Cadence1 ) )
        ihigh = int( numpy.round( (fhigh - self.fcoarse.Offset1) /
                                  self.fcoarse.Cadence1 ) )
        data = data[ ilow : ihigh + 1 ]
        fdata = self.fcoarse.data[ ilow : ihigh + 1 ]
        fscale = {'Offset1'  : fdata[0] ,
                  'Cadence1' : fdata[1] - fdata[0] }
        f             = Coarsable( fdata , **fscale )
        mlmlintegrand = Coarsable( data  , **fscale )
        return f , mlmlintegrand

    def getsummand( self , flow=1.4e-4 , fhigh=1e-1 , lmax=20 ):
        """Return the summand in the expression for the Fisher matrix
        Keyword parameter:
        flow -- lower frequency limit
        fhigh -- higher frequency limit
        Output:
        summand -- the summand
        """
        print 'Calculating GG...',
        if not hasattr( self , 'fcoarse' ) :
            self.getCoarseFreqs()
        if fhigh < flow :
            print "fhigh must be >= flow"
            raise ValueError
        elif flow < self.fcoarse.data[0] :
            print "flow must be >= %f" % self.fcoarse.data[0]
            raise ValueError
        elif fhigh > self.fcoarse.data[-1] :
            print "fhigh must be <= %f" % self.fcoarse.data[-1]
            raise ValueError
        
        glm = self.orf.getMultipleMoments( lmax=lmax )
        pII = self.psd['PowerSpectra'][ 'P'+2*self.pair[0] ]
        pJJ = self.psd['PowerSpectra'][ 'P'+2*self.pair[1] ]

        inputs = { 'Offset1'  : self.fcoarse.Offset1 ,
                   'Cadence1' : self.fcoarse.Cadence1 ,
                   'N1'       : len( self.fcoarse.data  ) }
        glm = glm.coarsegrain( **inputs )
        pII , pJJ = pII.coarsegrain( **inputs ) , pJJ.coarsegrain( **inputs )

        H = self.fcoarse.data**self.alpha

        ilow  = int( numpy.round( (flow - self.fcoarse.Offset1) /
                                  self.fcoarse.Cadence1 ) )
        ihigh = int( numpy.round( (fhigh - self.fcoarse.Offset1) /
                                  self.fcoarse.Cadence1 ) )

########## Option (1) Estimator with +ve and -ve frequencies ##############
#        H2oPP = H**2 / ( numpy.real( pII.data ) * numpy.real( pJJ.data ) )
#
#        lmax = self.orf.ntrunc
#        indxpn = LISAresponse.getMLvec( lmax , 'pn' )
#        Nml = len( indxpn )
#
#        Npml = int( 0.5 * ( Nml - ( lmax + 1 ) ) )
#        glmmdata = numpy.zeros( glm.data.shape , complex )
#        glmmdata[ : lmax+1 ] = glm.data[ : lmax+1 ]
#        glmmdata[ lmax+1 : lmax+1 + Npml ] = glm.data[ lmax+1 + Npml : ]
#        glmmdata[ lmax+1 + Npml : ] = glm.data[ lmax+1 : lmax+1 + Npml ]
#        minus1tom = numpy.zeros( ( len( indxpn ), ) )
#        for i,ml in enumerate( indxpn ):
#            minus1tom[ i ] = ( -1 )**ml[ 0 ]
#        minus1tomTglmmdata = numpy.transpose( numpy.transpose( glmmdata ) * minus1tom )
#        ############### If numpy.dot() works properly ####################
#        #pdata = numpy.dot( numpy.conj( glm.data ) * H2oPP , numpy.transpose( glm.data ) )
#        #ndata = numpy.dot( minus1tomTglmmdata * H2oPP , numpy.transpose( numpy.conj( minus1tomTglmmdata ) ) )
#        #data = pdata + ndata
#        ############### If numpy.dot() doesn't work properly ##############
#        pmatrix1 = numpy.conj( glm.data ) * H2oPP
#        pmatrix2 = numpy.transpose( glm.data )
#        nmatrix1 = minus1tomTglmmdata * H2oPP
#        nmatrix2 = numpy.transpose( numpy.conj( minus1tomTglmmdata ) )
#
#        data = numpy.zeros( ( Nml , Nml ) , dtype='complex' )
#        for ii in range( Nml ) :
#            for jj in range( Nml ) :
#                data[ ii , jj ] = numpy.sum( pmatrix1[ ii , ilow:ihigh+1 ] * pmatrix2[ ilow:ihigh+1 , jj ] ) + numpy.sum( nmatrix1[ ii , ilow:ihigh+1 ] * nmatrix2[ ilow:ihigh+1 , jj ] )
##############################################################################

################## Option (2) : loop over frequencies ########################
        H2oPP = H**2 / ( numpy.real( pII.data ) * numpy.real( pJJ.data ) )

        indxpn = LISAresponse.getMLvec( lmax , 'pn' )
        Nml = len( indxpn )

        Npml = int( .5 * ( Nml - ( lmax + 1 ) ) )
        glmmdata = numpy.zeros( glm.data.shape , complex )
        glmmdata[ : lmax+1 ] = glm.data[ : lmax+1 ]
        glmmdata[ lmax+1 : lmax+1 + Npml ] = glm.data[ lmax+1 + Npml : ]
        glmmdata[ lmax+1 + Npml : ] = glm.data[ lmax+1 : lmax+1 + Npml ]

        minus1tom = numpy.zeros( ( Nml, ) )
        for i,ml in enumerate( indxpn ) :
            minus1tom[ i ] = ( -1 )**ml[0]

        minus1tomTglmmdata = numpy.transpose( numpy.transpose( glmmdata ) * minus1tom )

        data = numpy.zeros( ( Nml , Nml ) , dtype='complex' )
        "~~ (l,m,l',m')s to check below"
        mlmls = [ ( (0,0),(0,0) ) , ( (-2,2),(1,5) ), ( (4,6),(-3,4) ) , ( (2,14),(-6,12) ) ]
        mlmlsm = [ ( ( -mlml[0][0] , mlml[0][1] ) , ( -mlml[1][0] , mlml[1][1] ) ) for mlml in mlmls ]
        inds  = [ ( indxpn.index( mlml[0] ) , indxpn.index( mlml[1] ) ) for mlml in mlmls ]
        indsm = [ ( indxpn.index( mlml[0] ) , indxpn.index( mlml[1] ) ) for mlml in mlmlsm ]
        Nf = ihigh - ilow + 1
        cheshs = [ numpy.zeros( (Nf,) , dtype='complex' ) for k in range( len( mlmls ) ) ]
        "~~"
        for ii in range( ilow , ihigh+1 ) :
            p = numpy.conj( glm.data[ : , ii ] )
            q = minus1tomTglmmdata[ : , ii ]
            pp = numpy.outer( p , numpy.conj(p) )
            qq = numpy.outer( q , numpy.conj(q) )
            ddata = H2oPP[ii] * ( pp + qq )
            data += ddata
            
            "~~ Get value for some (lm,l'm') at this frequency"
            for k in range( len( mlmls ) ) :
                cheshs[k][ ii - ilow ] = ddata[ inds[k] ]  
            "~~"
            
#            print "Constituent matrix 'pp(f)':"
#            if not IsHermitian( pp ) :
#                print "   " , "Warning: not Hermitian"
#            if not IsPositiveDefinite( pp ) :
#                print "   " , "Warning: not positive-definite."
#            if not IsUeqV( pp ) :
#                print "   " , "Warning: U != V in SVD"
#            print "   " , "condition number:" , ConditionNumber( pp )
#            print "   " , "smallest diagonals:" , SmallestDiagonals( pp )
#            U , s , Vh = pylab.svd( pp )
#            sinv = 1 / s
#            ppinv = numpy.dot( adjoint(Vh) , numpy.dot( numpy.diag(sinv) , adjoint(U) ) )
#            print "Its inverse:"
#            if not IsHermitian( ppinv ) :
#                print "   " , "Warning: not Hermitian"
#            if not IsPositiveDefinite( ppinv ) :
#                print "   " , "Warning: not positive-definite."
#            if not IsUeqV( ppinv ) :
#                print "   " , "Warning: U != V in SVD"
#            print "   " , "condition number:" , ConditionNumber( ppinv )
#            print "   " , "smallest diagonals:" , SmallestDiagonals( ppinv )
#            if not InvertWell( ppinv ) :
#                print "   " , "Warning: doesn't invert well."
##########################################################################
        "~~ Write H2oPP, cheshs and glms to disk"
        cheshdirs = [ 'cheshs/chesh1/' , 'cheshs/chesh2/' , 'cheshs/chesh3/' , 'cheshs/chesh4/' ]
        for cheshdir in cheshdirs :
            if cheshdir not in glob.glob( cheshdir ) :
                os.system( 'mkdir -p %s' % cheshdir )
        
        f_chesh = self.fcoarse.data[ ilow:ihigh+1 ]
        cheshdicts = [ { 'f' : f_chesh ,
                         'H2oPP' : H2oPP[ ilow:ihigh+1 ] ,
                         'chesh' : cheshs[k] ,
                         'ml1' : mlmls[k][0] , 'ml2' : mlmls[k][1] ,
                         'glm1' : glm.data[ inds[k][0] , ilow:ihigh+1 ] ,
                         'glm2' : glm.data[ inds[k][1] , ilow:ihigh+1 ] ,
                         'glm1m' : glm.data[ indsm[k][0] , ilow:ihigh+1 ] ,
                         'glm2m' : glm.data[ indsm[k][1] , ilow:ihigh+1 ]
                         }
                       for k in range( len(mlmls) ) ]
        
        for k in range( len( mlmls ) ) :
            file = open( cheshdirs[k]+'chesh-s%03d.pkl' % self.day , 'wb' );
            cPickle.dump( cheshdicts[k] , file , -1 )
            file.close()
        "~~"
        summand = Coarsable( data )

        print 'done'
        return summand



def Adjoint( matrix ) :
    return numpy.transpose( numpy.conj( matrix ) )


def IsPositiveDefinite( A , smallest=0 ) :
    """
    INPUT:
    A --- matrix ( numpy array )
    smallest --- print  the smallest 'smallest' Evalues
    OUTPUT:
    Boolean --- True or False
    """
    d , P = pylab.eig( A )
    posdef = True
    for ii in range( d.shape[0] ) :
        if numpy.real( d[ ii ] ) <= 0 :
            posdef = False
            break
    if smallest :
        print numpy.sort( d )[ :smallest ]
    return posdef


def IsUeqV( A , abscomp=False ) :
    """
    INPUT:
    A --- matrix (numpy array)
    abscomp --- Boolean, if True returns (abs of mean) and (max abs of diff)
    OUTPUT:
    Boolean --- True or False
    """
    U , s , Vh = pylab.svd( A )
    if abscomp :
        max_abs_diff = numpy.max( numpy.abs( U - Adjoint(Vh) ) )
        abs_mean = numpy.mean( numpy.abs( U ) )
        print "max. abs(diff):" , max_abs_diff
        print "mean of abs(matrix):" , abs_mean
    return numpy.allclose( U , Adjoint(Vh) )


def DiagonalsAllPositive( A , smallest=0 ) :
    """
    INPUT:
    A --- matrix (numpy array)
    OUTPUT:
    Boolean --- True or False
    """
    d = numpy.diag( A )
    allpos = True
    for ii in range( d.shape[0] ) :
        if numpy.real( d[ ii ] ) <= 0 :
            allpos = False
            break
    if smallest :
        print numpy.sort( d )[ :smallest ]
    return allpos


def ConditionNumber( A ) :
    """
    INPUT:
    A --- matrix ( numpy array )
    OUTPUT:
    cn --- condition number of A
    """
    U , s , Vh = pylab.svd( A )
    cn = numpy.max(s) / numpy.min(s)
    return cn


def InvertWell( A , abscomp=False ) :
    """
    INPUT:
    A --- matrix ( numpy array )
    abscomp --- Boolean, if True returns (abs of mean) and (max abs of diff)
    OUTPUT:
    Boolean --- True or False
    """
    U , s , Vh = pylab.svd( A )
    sinv = 1 / s
    Ainv = numpy.dot( Adjoint(Vh) , numpy.dot( numpy.diag(sinv) , Adjoint(U) ) )
    AinvA = numpy.dot( Ainv , A )
    I = numpy.identity( A.shape[0] )
    if abscomp :
        max_abs_diff = numpy.max( numpy.abs( AinvA - I ) )
        abs_mean = numpy.mean( numpy.abs( AinvA ) )
        print "max. abs(diff):" , max_abs_diff
        print "mean of abs(matrix):" , abs_mean
    return numpy.allclose( AinvA , I )

    
def IsHermitian( matrix , abscomp = False ) :
    """
    INPUT:
    matrix --- numpy array
    abscomp --- Boolean, if True returns (abs of mean) and (max abs of diff)
    OUTPUT:
    Boolean --- True or False
    """
    if abscomp :
        max_abs_diff = numpy.max( numpy.abs( matrix - Adjoint(matrix) ) )
        abs_mean = numpy.mean( numpy.abs( matrix ) )
        print "max. abs(diff):" , max_abs_diff
        print "mean of abs(matrix):" , abs_mean
    return numpy.allclose( Adjoint(matrix) , matrix )


def findTSfiles(sourcename,tsdir,day):
    tsfileprefix = sourcename + '_obs-s%03d' % day
    tsxmlname = tsfileprefix + '.xml'
    tsbinname = tsfileprefix + '-0.bin'
    workdir = os.getcwd()
    os.chdir(tsdir)
    for name in [tsxmlname,tsbinname]:
        if name in glob.glob(tsfileprefix + '*'):
            available = True
        else:
            available = False
            print ( name + ' not found in ' + os.getcwd() )
            xmlpath = []
            break
    else:
        xmlpath = [tsdir + tsxmlname]
    os.chdir(workdir)
    return available, xmlpath

def findpklTSfiles( sourcename , tsdir , day ):
    tsnameprefix = sourcename + '_obs-s%03d' % day
    tsname = tsnameprefix + '.pkl'
    workdir = os.getcwd()
    os.chdir( tsdir )
    if tsname in glob.glob( tsnameprefix + '*' ):
        available = True
        tspath = [ tsdir + tsname ]
    else:
        available = False
        print ( tsname + ' not found in ' + os.getcwd() )
        tspath = []

    os.chdir( workdir )
    return available , tspath
    

def findPSDfiles(sourcename,psddir,day):
    psdnameprefix = sourcename + '_obs-s%03d-psd' % day
    psdname = psdnameprefix + '.pkl'
    workdir = os.getcwd()
    os.chdir( psddir )
    for name in [psdname]:
        if name in glob.glob( psdnameprefix + '*' ):
            available = True
        else:
            available = False
            print ( name + ' not found in ' + os.getcwd() )
            psdpath = []
            break
    else:
        psdpath = [ psddir + psdname ]
    os.chdir( workdir )
    return available , psdpath
    

def findORFfiles(orfdir,day):
    orfprefix = 'orf-d%03d' % day
    orfname = orfprefix + '.pkl'
    workdir = os.getcwd()
    os.chdir( orfdir )
    for name in [ orfname ] :
        if name in glob.glob( orfprefix + '*' ):
            available = True
        else:
            available = False
            print name , 'not found in' , os.getcwd()
            orfpath = []
            break
    else:
        orfpath = [ orfdir + name for name in [ orfname ] ]
    os.chdir( workdir )
    return available , orfpath


def findXXfiles( sourcename , XXdir , day ):
    XXnameprefix = sourcename + '_obs-s%03d-XX' % day
    XXname = XXnameprefix + '.pkl'
    workdir = os.getcwd()
    os.chdir( XXdir )
    for name in [ XXname ]:
        if name in glob.glob( XXnameprefix + '*' ):
            available = True
        else:
            available = False
            print ( name + ' not found in ' + os.getcwd() )
            XXpath = []
            break
    else:
        XXpath = [ XXdir + XXname ]
    os.chdir( workdir )
    return available , XXpath

def findGGfiles( sourcename , GGdir , day ):
    GGnameprefix = sourcename + '_obs-s%03d-GG' % day
    GGname = GGnameprefix + '.pkl'
    workdir = os.getcwd()
    os.chdir( GGdir )
    for name in [ GGname ]:
        if name in glob.glob( GGnameprefix + '*' ):
            available = True
        else:
            available = False
            print ( name + ' not found in ' + os.getcwd() )
            GGpath = []
            break
    else:
        GGpath = [ GGdir + GGname ]
    os.chdir( workdir )
    return available , GGpath
                                                                                                        

def tsXML_to_psdPKL(tspath,psddir):

    tsdict = tsXML_to_tsDICT( tspath )
    ts = TimeSeries( tsdict ) 
    
    segduration = 2*60.0**2
    nfft = int( segduration / ts.t.Cadence1 )

    f , PAA = ts.psd( 'A' , nfft , coarsable=True )
    f , PEE = ts.psd( 'E' , nfft , coarsable=True )
    f , PTT = ts.psd( 'T' , nfft , coarsable=True )

    psddic = { 'PowerSpectra' : 
              {'f': f , 'PAA': PAA , 'PEE': PEE , 'PTT': PTT} }

    tsname = os.path.basename(tspath)
    psdname = re.sub( '\.xml$' , '' , tsname ) + '-psd.pkl'
    psdpath = psddir + psdname
    file = open( psdpath , 'wb' )
    cPickle.dump( psddic , file , -1 )
    file.close()
    return

def tsDICT_to_psdPKL( tsdict , psdpath ):

    ts = TimeSeries( tsdict )

    segduration = 2 * 60.**2
    nfft = int( segduration / ts.t.Cadence1 )

    f , PAA = ts.psd( 'A' , nfft , coarsable=True )
    f , PEE = ts.psd( 'E' , nfft , coarsable=True )
    f , PTT = ts.psd( 'T' , nfft , coarsable=True )

    psddict = {}
    psddict['PowerSpectra'] = { 'f': f ,
                                'PAA': PAA , 'PEE': PEE , 'PTT': PTT }

    file = open( psdpath , 'wb' )
    cPickle.dump( psddict , file , -1 )
    file.close()
    return


def avg_psdPKL_avg(day,lpsdpath,rpsdpath,psddir):
    file = open( lpsdpath , 'rb' )
    lpsd = cPickle.load( file )
    file.close()
    file = open( rpsdpath , 'rb' )
    rpsd = cPickle.load( file )
    file.close()
    
    PAA=(lpsd['PowerSpectra']['PAA'].data + rpsd['PowerSpectra']['PAA'].data)/2
    PEE=(lpsd['PowerSpectra']['PEE'].data + rpsd['PowerSpectra']['PEE'].data)/2
    PTT=(lpsd['PowerSpectra']['PTT'].data + rpsd['PowerSpectra']['PTT'].data)/2

    f = lpsd['PowerSpectra']['f']
    PAA = Coarsable( PAA , Offset1=f.Offset1 , Cadence1=f.Cadence1 )
    PEE = Coarsable( PEE , Offset1=f.Offset1 , Cadence1=f.Cadence1 )
    PTT = Coarsable( PTT , Offset1=f.Offset1 , Cadence1=f.Cadence1 )

    psd = {}
#    psd['LISAData'] = lpsd['LISAData']
#    psd['SourceData'] = lpsd['SourceData']
    psd['PowerSpectra'] = { 'f':f , 'PAA':PAA , 'PEE':PEE , 'PTT':PTT }

    lpsdname = os.path.basename( lpsdpath )
    psdname = re.split( '-' , lpsdname , 1 )[0] + '-s%03d-psd.pkl' % day
    psdpath = psddir + psdname
    file = open( psdpath , 'wb' )
    cPickle.dump( psd , file , -1 )
    file.close()
    return



from mpl_toolkits.basemap import Basemap , addcyclic
class SkyMap(object):
    def __init__( self , xlmpath , lmax=15 ):
        file = open( xlmpath , 'rb' )
        xlmdic = cPickle.load( file )
        file.close()
        self.maptype         = xlmdic['maptype']
#        self.psdData         = xlmdic['psdData']
        self.orfData         = xlmdic['orfData']
        self.GWSpectralSlope = xlmdic['GWSpectralSlope']
        self.ntrunc          = xlmdic['ntrunc']
        self.xlm = get_lmax_subset_from( xlmdic['MultipleMoments'].data , lmax )
        self.gotpix = False
        xlmname = os.path.basename( xlmpath )
        self.outfileprefix = '-'.join( re.split( '-' , xlmname )[:-1] )

        #Quick handling of available 'singular values'
        if xlmdic.has_key( 'Singularvalue' ):
            self.singular_value = xlmdic[ 'Singularvalue' ]

        return

    def getMultipleMoments( self , msign = 'pn' ):
        if msign == 'p' :
            indxp = LISAresponse.getMLvec( self.ntrunc , 'p' )
            xlm = self.xlm[ : len(indxp) ]
        elif msign == 'pn' :
            xlm = self.xlm
        else:
            print "Keyword parameter 'msign' must either be 'p' for only positive m's or 'pn' for all m's."
            raise
        return xlm

    def getRealMultipleMoments( self , msign = 'pn' ):
        indxpn = LISAresponse.getMLvec( self.ntrunc , 'pn' )
        plm = numpy.zeros( numpy.shape( self.xlm ) , complex )
        for i,ml in enumerate( indxpn ) :
            m , l = ml[0] , ml[1]
            k = indxpn.index( ( -m , l ) )
            plm[i] = ( self.xlm[i] + (-1)**m*numpy.conj(self.xlm[k]) ) / 2
        if msign == 'p' : 
            indxp = LISAresponse.getMLvec( self.ntrunc , 'p' )
            plm = plm[ : len( indxp ) ]
        elif msign == 'pn' :
            plm = plm
        else:
            print "Keyword parameter 'msign' must either be 'p' for only positive m's or 'pn' for all m's."
            raise
        return plm

    def getImagMultipleMoments( self , msign = 'pn' ):
        indxpn = LISAresponse.getMLvec( self.ntrunc , 'pn' )
        qlm = numpy.zeros( numpy.shape( self.xlm ) , dtype = complex )
        for i,ml in enumerate( indxpn ):
            m , l = ml[0] , ml[1]
            k = indxpn.index( ( -m , l ) )
            qlm[i] = ( self.xlm[i] - (-1)**m*numpy.conj( self.xlm[k] ) ) / 2j
        if msign == 'p' : 
            indxp = LISAresponse.getMLvec( self.ntrunc , 'p' )
            qlm = qlm[ : len(indxp) ]
        elif msign == 'pn' : 
            qlm = qlm
        else:
            print "Keyword parameter 'msign' must either be 'p' for only positive m's or 'pn' for all m's."
            raise
        return qlm

    def getPixels( self ):
        nlat = self.ntrunc + 50
        nlon = 2*( nlat - 1 )
        self.sky = LISAresponse.mySpharmt( nlon , nlat )
        if self.maptype in [ 'X' , 'XX' ] :
            norm = 1e-45
        elif self.maptype in [ 'P' ] :
            norm = 1e45
        elif self.maptype in [ 'K' , None ] :
            norm = 1

        indxp = LISAresponse.getMLvec( self.ntrunc , 'p' )
        to_pyspharm = numpy.zeros( ( len(indxp), ) )
        for i , ml in enumerate( indxp ):
            to_pyspharm[i] = (-1)**ml[0] / numpy.sqrt( 2*numpy.pi )
        plm = to_pyspharm * self.getRealMultipleMoments( msign='p' ) * norm
        qlm = to_pyspharm * self.getImagMultipleMoments( msign='p' ) * norm

        self.pixp = self.sky.spectogrd( plm )
        self.pixq = self.sky.spectogrd( qlm )
        self.gotpix = True
        return self.pixp , self.pixq


    def plotMap( self , figdir='' ):
        if not self.gotpix:
            self.getPixels()

        lats = self.sky.lats
        lons = self.sky.lons
        pixpw , lonsw = addcyclic( self.pixp , lons )
        pixqw , lonsw = addcyclic( self.pixq , lons )
        meshlon , meshlat = pylab.meshgrid( lonsw , lats )
        
        projection = Basemap( projection='moll' , lon_0=180 , resolution='c' )
        projection.drawmapboundary()
        x , y = projection( meshlon , meshlat )

        fig = pyplot.figure( figsize = (8,8) )
        ax = fig.add_axes( [ 0.05 , 0.15 , 0.8 , 0.8] )
        projection.contourf( x , y , pixpw , 30 )
        projection.drawmeridians( numpy.arange(0,360,30) )
        projection.drawparallels( numpy.arange(-90,90,30) , labels=[1,0,0,0] )
        ####Quick labelling of available singular value
        if hasattr( self , 'singular_value' ):
            pyplot.title( '%e' % self.singular_value )
        ####
        pos = ax.get_position()
        l , b , w , h = pos.bounds
        cax = pyplot.axes( [ l+w+0.03 , b , 0.04 , h ] )
        pyplot.colorbar( cax=cax , orientation='vertical' )
        figname = self.outfileprefix + '-' + str( self.maptype ) + 'real.png'
        pylab.savefig( figdir + figname )

        fig = pyplot.figure( figsize = (8,8) )
        ax = fig.add_axes( [ 0.05 , 0.15 , 0.8 , 0.8 ] )
        projection.contourf( x , y , pixqw , 30 )
        projection.drawmeridians( numpy.arange(0,360,30) )
        projection.drawparallels( numpy.arange(-90,90,30) , labels=[1,0,0,0])
        ####Quick labelling of available singular value
        if hasattr( self , 'singular_value' ):
            pyplot.title( '%e' % self.singular_value )
        ####
        pos = ax.get_position()
        l , b , w , h = pos.bounds
        cax = pyplot.axes( [ l+w+0.03 , b , 0.04 , h ] )
        pyplot.colorbar( cax=cax , orientation='vertical' )
        figname = self.outfileprefix + '-' + str( self.maptype ) + 'imag.png'
        pylab.savefig(  figdir + figname  )
        return



def IsMinusOneP( matrix , lmax ) :
    indxpn = LISAresponse.getMLvec( lmax , 'pn' )
    Nml = len( indxpn )
    Npml = int( 0.5 * ( Nml - ( lmax + 1 ) ) )

    minus1tomPm = numpy.zeros( ( Nml , Nml ) )
    for i1,ml1 in enumerate( indxpn ) :
        for i2 , ml2 in enumerate( indxpn ) :
            minus1tomPm[ i1 , i2 ] = ( -1 )**( ml1[0] + ml2[0] )
    datamm = numpy.zeros( matrix.shape , dtype = matrix.dtype )
    datamm[ :lmax+1 , :lmax+1 ]            = matrix[ :lmax+1 , :lmax+1 ]
    datamm[ :lmax+1 , lmax+1:lmax+1+Npml ] = matrix[ :lmax+1 , lmax+1+Npml: ]
    datamm[ :lmax+1 , lmax+1+Npml: ]       = matrix[ :lmax+1 , lmax+1:lmax+1+Npml ]
    datamm[ lmax+1:lmax+1+Npml , :lmax+1 ]            = matrix[ lmax+1+Npml: , :lmax+1 ]
    datamm[ lmax+1:lmax+1+Npml , lmax+1:lmax+1+Npml ] = matrix[ lmax+1+Npml: , lmax+1+Npml: ]
    datamm[ lmax+1:lmax+1+Npml , lmax+1+Npml: ]       = matrix[ lmax+1+Npml: , lmax+1:lmax+1+Npml ]
    datamm[ lmax+1+Npml: , :lmax+1 ]            = matrix[ lmax+1:lmax+1+Npml , :lmax+1 ]
    datamm[ lmax+1+Npml: , lmax+1:lmax+1+Npml ] = matrix[ lmax+1:lmax+1+Npml , lmax+1+Npml: ]
    datamm[ lmax+1+Npml: , lmax+1+Npml: ]       = matrix[ lmax+1:lmax+1+Npml , lmax+1:lmax+1+Npml ]
    minus1tomPmTdatamm = minus1tomPm * datamm
    return numpy.allclose( minus1tomPmTdatamm , numpy.conj( matrix ) )

    



class FisherMatrix(object):
    def __init__( self , fishpath , lmax=15 ):
        file = open( fishpath , 'rb' )
        fishdic = cPickle.load( file )
        file.close()
        self.fishtype        = fishdic['fishtype']
        self.orfData         = fishdic['orfData']
        self.GWSpectralSlope = fishdic['GWSpectralSlope']
        self.ntrunc          = lmax 
        self.fish = get_lmax_subset_from( fishdic['FisherMatrix'].data , lmax )

        if not IsHermitian( self.fish ) :
            print "Warning: Fisher matrix NOT Hermitian."

        if not IsMinusOneP( self.fish , self.ntrunc ) :
            print "Warning: Fisher matrix does not satisfy the minus 1 to the m+m' identity."

        self.decomposed  = False
        self.regularised = False
        return


    def svd(self):
        self.U, self.s, self.Vh = pylab.svd(self.fish)
        self.decomposed = True
        return self.U , self.s , self.Vh


    def regularise(self,regMethod=2,regCutoff=2./3):
        if not self.decomposed :
            self.svd()
        if regMethod in [ 1 , 11 ]:
            ind = pylab.find( self.s/self.s[0] >= regCutoff )
            if len(ind) > 0 :
                ii = ind[-1]
            else:
                ii = 0
        elif regMethod in [ 2 , 22 ]:
            ii = int( numpy.round( len(self.s)*regCutoff ) ) - 1
        elif regMethod in [ 3 , 33 ]:
            ii = int( numpy.round( regCutoff ) ) - 1
        else:
            print 'Unknown regularisation method.'
            raise Exception
        if ii < 1 :
            ii = 0
            if ii > len(self.s)-1 :
                ii = len(self.s)-1 #-1
        self.regs = numpy.zeros(self.s.shape,self.s.dtype)
        if regMethod in [ 1 , 2 , 3 ]:
            self.regs[:ii+1] , self.regs[ii+1:] = self.s[:ii+1] , self.s[ii]
        elif regMethod in [ 11 , 22 , 33 ]:
            self.regs[:ii+1] , self.regs[ii+1:] = self.s[:ii+1] , numpy.inf
        else:
            print 'Unknown regularisation method.'
        self.regularised = True
        self.N_keptSV = ii + 1
        return self.regs


    def invert(self):
        if not self.decomposed:
            self.svd()
        self.sinv = 1. / self.s
        Sinv = numpy.diag( self.sinv )
        Uh = numpy.transpose(numpy.conj(self.U))
        V  = numpy.transpose(numpy.conj(self.Vh))
        fishinv = numpy.dot( V , numpy.dot( Sinv , Uh ) )
        for ii in range( fishinv.shape[0] ) :
            fishinv[ ii , ii ] = numpy.real( fishinv[ ii , ii ] )
        self.fishinv = fishinv

        if not IsHermitian( self.fishinv ) :
            print "Warning: Fisher matrix's inverse NOT Hermitian."
        if not IsMinusOneP( self.fishinv , self.ntrunc ) :
            print "Warning: Fisher matrix's inverse does not satisfy the minus 1 to the m+m' identity."
        return self.fishinv

    def reginvert(self):
        if not self.regularised :
            self.regularise()
        self.regsinv = 1 / self.regs
        Sinv = numpy.diag( self.regsinv )
        Uh = numpy.transpose(numpy.conj(self.U))
        V  = numpy.transpose(numpy.conj(self.Vh))
        regfishinv = numpy.dot( V , numpy.dot( Sinv , Uh ) )

        for ii in range( regfishinv.shape[0] ) :
            regfishinv[ ii , ii ] = numpy.real( regfishinv[ ii , ii ] )
        self.regfishinv = regfishinv

        if not IsHermitian( self.regfishinv ) :
            print "Warning: Fisher matrix's inverse NOT Hermitian."
        if not IsMinusOneP( self.regfishinv , self.ntrunc ) :
            print "Warning: Fisher matrix's inverse does not satisfy the minus 1 to the m+m' identity."
        return self.regfishinv


    def getpCovar(self):
        if not self.regularised :
            self.regularise()
        Cov = numpy.diag( self.s / self.regs**2 )
        Uh = numpy.transpose(numpy.conj(self.U))
        V  = numpy.transpose(numpy.conj(self.Vh))
        self.pCovar = numpy.dot(V,numpy.dot(Cov,Uh))
        return self.pCovar
                                                            


def Deconvolve( skymap , fishermatrix , regMethod=2 , regCutoff=2./3 ):
    if skymap.orfData['Antenna'] != fishermatrix.orfData['Antenna'] :
        print "Skymap and Fisher matrix are from different detector pair"
        raise
    if skymap.ntrunc != fishermatrix.ntrunc :
        print "lmax of skymap and Fisher matrix are different"
        raise Exception , "lmax of skymap and Fisher matrix are different"
    if not regMethod:
        fishinv = fishermatrix.invert()
    else:
        fishermatrix.regularise( regMethod , regCutoff )
        fishinv = fishermatrix.reginvert()
    cleanmap = numpy.dot( fishinv , skymap.xlm )
    return cleanmap



def get_MLvec_indices( lmax0 , lmax ) :
    """
    INPUT:
    lmax0 --- original lmax
    lmax  --- new lmax(should be smaller than lmax0)
    OUTPUT:
    indices --- indices of (m,l)'s from mlvec in mlvec0
    """
    if lmax > lmax0 :
        raise ValueError , "New lmax must be smaller or equal to old lmax."
    
    mlvec0 = LISAresponse.getMLvec( lmax0 )
    mlvec = LISAresponse.getMLvec( lmax )
    
    indices = []
    for ml in mlvec :
        indices += [ mlvec0.index( ml ) ]
    return indices

def get_MLvec_indices_from_l( lmax0 , l ) :
    """
    INPUT:
    lmax0 --- original lmax
    l     --- l to restric lmax to
    OUTPUT:
    indices --- indices of (m,l) in mlvec0
    """
    if l > lmax0 :
        raise ValuesError , "l must be smaller or equal to lmax0."

    mlvec0 = LISAresponse.getMLvec( lmax0 )

    indices = []
    for i,ml in enumerate( mlvec0 ) :
        if ml[-1] == l :
            indices += [i]
        else :
            continue
    return indices 
            
    
def get_lmax_subset_from( matrix , lmax , axis='All' ) :

    if axis == 0 :
        A = numpy.copy( matrix )
        if numpy.mod( numpy.sqrt( A.shape[0] ) - 1 , 1 ) != 0 :
            raise Exception , "Side must have length corresponding to an integer lmax."
        else:
            lmax0 = int( numpy.sqrt( A.shape[0] ) - 1 )
        
        indices = get_MLvec_indices( lmax0 , lmax )
        new_A = A[ indices ]
        return new_A

    elif axis == 1 :
        A = numpy.transpose( numpy.copy( matrix ) )
        new_A = numpy.transpose( get_lmax_subset_from( A , lmax , axis=0 ) )
        return new_A

    elif axis == 'All' :
        A = get_lmax_subset_from( numpy.copy( matrix ) , lmax , axis=0 )
        new_A = get_lmax_subset_from( A , lmax , axis=1 )
        return new_A
    
    else:
        raise ValueError , "axis must either be 0, 1, or 'All'"


    
def get_l_subset_from( matrix , l , axis='All' ) :
    
    if axis == 0 :
        A = numpy.copy( matrix )
        if numpy.mod( numpy.sqrt( A.shape[0] ) - 1 , 1 ) != 0 :
            raise Exception , "Side must have length corresponding to an integer lmax."
        else:
            lmax0 = int( numpy.sqrt( A.shape[0] ) - 1 )

        indices = get_MLvec_indices_from_l( lmax0 , l )
        new_A = A[ indices ]
        return new_A

    elif axis == 1 :
        A = numpy.transpose( numpy.copy( matrix ) )
        new_A = numpy.transpose( get_l_subset_from( A , l , axis=0 ) )
        return new_A

    elif axis == 'All' :
        A = get_l_subset_from( numpy.copy( matrix ) , l , axis=0 )
        new_A = get_l_subset_from( A , l , axis=1 )
        return new_A
    
    else:
        raise ValueError , "axis must either be 0, 1, or 'All'"



def covarSpH2Pixel( covar , nlat , nlon ) :
    """
    INPUT:
    covar --- covariance matrix in SpH basis
    nlat --- number of latitudes
    nlon --- number of longitudes
    OUTPUT:
    covarP --- covariance matrix in pixel basis
    lats --- latitudes
    lons --- longitudes
    """
    lmax = int( numpy.sqrt( covar.shape[0] ) - 1 )
    U , lats , lons = U4.checkPIXELCONVERSION( lmax , nlat , nlon )
    Uh = numpy.transpose( numpy.conj( U ) )
    covarP = numpy.dot( U , numpy.dot( covar , Uh ) )
    return covarP , lats , lons


def getSigmaMap( covar , nlat , nlon ) :
    """
    INPUT:
    covar --- covariance matrix in SpH basis
    nlat --- number of latitudes
    nlon --- number of longitudes
    OUTPUT:
    sigmaP --- standard deviation skyap in pixel basis
    lats --- latitudes
    lons --- longitudes
    """
    covarP , lats , lons = covarSpH2Pixel( covar , nlat , nlon )
    sigmaP = numpy.sqrt( numpy.real( numpy.diag( covarP ) ) )
    return sigmaP , lats , lons



def order_ML_to_LM( ml_matrix ) :
    """
    Re-orders a matrix whose rows are ordered as the (m,l)s in getMLvec(),
    in the order l**2 + l + m + 1 
    INPUT:
    ml_matrix --- numpy array whose rows are to be re-ordered
    OUTPUT:
    lm_matrix --- numpy array whose rows have been re-ordered
    """
    Nml = ml_matrix.shape[0]
    lmax = numpy.sqrt( Nml ) - 1
    if ( lmax % 1 ) != 0 :
        raise Exception , 'Number of rows does not correspond to an integer lmax.'
    else:
        lmax = int( lmax )
        
    indxpn = LISAresponse.getMLvec( lmax )

    ml_matrixl = list( ml_matrix )
    lm_matrixl = []
    for l in range( lmax+1 ) :
        for m in range( -l , l+1 ) :
            lm_matrixl += [ ml_matrixl[ indxpn.index( (m,l) ) ] ]
    lm_matrix = numpy.array( lm_matrixl )
    return lm_matrix


def get_Angular_PSD( mla ) :
    """
    Returns angular PSD from an array whose elements are arranged according to (m,l) from LISAresponse.getMLvec()
    INPUT:
    mla --- multipole moments aranged by LISAresonse.getMLvec( msign='pn' ). (Contains m<0)
    OUTPUT:
    C --- angular power spectral density
    """
    Nml = mla.shape[0]
    lmax = numpy.sqrt( Nml ) - 1
    if ( lmax % 1 ) != 0 :
        raise Exception , 'Number of rows does not correspond to an integer lmax.'
    else:
        lmax = int( lmax )

    indxpn = LISAresponse.getMLvec( lmax )

    mlal = list( mla )

    C = numpy.zeros( ( lmax+1 , ) )
    for l in range( lmax+1 ) :
        for m in range( -l , l+1 ) :
            mlv = mlal[ indxpn.index( (m,l) ) ]
            C[l] += numpy.abs( mlv )**2
        C[l] = C[l] / ( 2*l + 1 ) 
    return C
