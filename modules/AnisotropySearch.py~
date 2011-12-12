import numpy as np
import sys
import glob
import os
import re
import cPickle as cpkl
import cmath
import spharm 
import scipy as sp
import scipy.linalg as sp_linalg
import lisaxml







def getMLvec(ntrunc,m='pn'):
    """
    Returns a list of (m,l) = (order,degree) up to l = ntrunc
    Input:
    ntrunc -- maximum l
    m -- 'p' for 0 < m < l, or 'pn' for -l < m < l
    Output:
    A list where (m,l) are arranged as, letting L = ntrunc
    m | 0 0 ... 0 1 1 ... 1 ... L-1 L-1 L
    -------------------------------------  for 0 < m < l , 0 < l < L , or
    l | 0 1 ... L 1 2 ... L ... L-1   L L

    m | 0 0 ... 0 1 1 ... 1 ... L-1 L-1 L -1 -1 ... -1 ... -L+1 -L+1 -L 
    --------------------------------------------------------------------  for -l < m < l , 0 < l < L  .
    l | 0 1 ... L 1 2 ... L ... L-1   L L  1  2 ...  L ...  L-1  L-1  L  

    NOTE: The getMLvec function from PYSPHARM is used here basically.
    """
    if ntrunc==0:
        M = [0]
        L = [0]
    elif m=='pn':    
        m,l = spharm.getspecindx(ntrunc)
        M = list(m) + list(-m[ntrunc+1:])
        L = list(l) + list(l[ntrunc+1:])
    elif m=='p':
        M,L = spharm.getspecindx(ntrunc)
    return zip(M,L)


class mySpharmt( spharm.Spharmt ) :
    """
    Same as a spharm.Spharmt instance, but with additional attributes:
    lats --- array of latitudes on the sphere in degrees
    lons --- array of longitudes on the sphere in degrees
    """
    def __init__(self,nlon,nlat,rsphere=1e2,gridtype='regular',legfunc='stored'):
        spharm.Spharmt.__init__(self,nlon,nlat,rsphere,gridtype,legfunc)
        if self.gridtype=='gaussian':
            raise Exception, "Can't deal with Gaussian grid at the moment."
        else:
            deltalat = 180./(nlat-1)
            self.lats = 90 - deltalat*np.arange(nlat)
        deltalon = 360./nlon
        self.lons = deltalon*np.arange(nlon)
        return
    

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
    y = np.zeros(Ny,complex)
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
        data = np.array( [ data ] )
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
    i = np.arange(Ny)
    jlow = np.intp(
        np.ceil(
        np.round( ( flowy + (i-0.5)*deltafy  -  flowx + 0.5*deltafx ) / deltafx , 10 )
        ) - 1
        )
    jhigh = np.intp(
        np.floor(
        np.round( ( flowy + (i+0.5)*deltafy - flowx + 0.5*deltafx  ) / deltafx , 10 )
        )
        )
    index1 = jlow[0]
    index2 = jhigh[-1]
    fraclow  = ( flowx+(jlow+0.5)*deltafx - flowy-(i-0.5)*deltafy  ) / deltafx
    frachigh = ( flowy+(i+0.5)*deltafy - flowx-(jhigh-0.5)*deltafx ) / deltafx
    frac1 = fraclow[0]
    frac2 = frachigh[-1]
    jtemp = jlow + 1

    coarsedata = np.zeros( (data.shape[0],Ny) , dtype = data.dtype )
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

        coarsedata[lm,:] = np.array( list(ya) + [yb] )
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
#                    print "Trying to coarsegrain side 1..."
                    coarsedata = coarsegrain(self.data,self.Offset1,self.Cadence1,Offset1,Cadence1,N1)
#                    print "Side 1 coarsegrained."
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
                    data       = np.transpose(self.data)
                    coarsedata = coarsegrain(     data,self.Offset0,self.Cadence0,Offset0,Cadence0,N0)
                    print "Side 0 coarsegrained."
                    coarsedata = Coarsable( np.transpose(coarsedata) , 
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
        Nlow = np.ceil( (highest_led - (F.data[0]-F.Cadence1/2.)) / F.Cadence1 ) 

    if F.data[-1]+F.Cadence1/2. <= lowest_hed :
        data = F.data[ Nlow : ]
    else:
        Nhigh = np.floor( (lowest_hed - (F.data[0]+F.Cadence1/2.)) / F.Cadence1 ) 
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
    tsdict['s1'] = Coarsable( tdiobs.A , **ts_scale )
    tsdict['s2'] = Coarsable( tdiobs.E , **ts_scale )
    tsdict['s3'] = Coarsable( tdiobs.T , **ts_scale )
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
        tsdict -- dictionary containing 't','1','2','3'

        OUTPUT:
        A TimeSeries object containing time-series info

        """
        self.t = tsdict['t']
        self.A = tsdict['1']
        self.E = tsdict['2']
        self.T = tsdict['3']
        return
    def scale_by( self , factor ) :
        """
        scales the time-series by a constant factor
        INPUT:
        factor --- factor by which the time-seris is to be scaled
        """
        self.A.data *= factor
        self.E.data *= factor
        self.T.data *= factor
        return
#20110205: self.psd() and self.csd() removed because matplotlib/pylab cannot be imported on the computing nodes of the cluster.
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
        z = np.zeros( fftlength - len(x.data) )
        xbar = sp.fft( np.array( list(x.data*window) + list(z) ) )
        winfftx = xbar[ 0:np.floor( fftlength/2+1 ) ]*x.Cadence1
        df = 1. / (x.Cadence1*fftlength)
        f =  df*np.arange( len(winfftx) )
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
        norm = window.shape[0] / np.sum( window**2 ) 
        Cxy = np.conj(winfftx)*winffty * 2 / dataduration * norm
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
        orf  = cpkl.load( file )
        file.close()
#        self.LISAData = orf['LISAData']
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
    def getMultipleMoments(self,msign,lmax):
        """Return multiple moments of orf

        Keyword parameters:
        msign -- 'p'  -- return only positive m's multiple moments
              -- 'pn' -- return all multiple moments

        Output:
        g -- Coarsable object containing orf's multiple moments

        """
        sqrt2pi = np.sqrt( 2*np.pi )
        if msign == 'p' :
            print "Warning: input lmax ineffective for msign='p'"
            p , q = np.conj( self.real.data ) , np.conj( self.imag.data )
            indxp = getMLvec( self.ntrunc , 'p' )
            data = np.zeros( p.shape , dtype='complex' )
            for i,ml in enumerate( indxp ):
                data[i] = (-1)**ml[0] * sqrt2pi * ( p[i] + 1j*q[i] )
            g = Coarsable( data , **self.fscale )
        elif msign == 'pn' :
            p , q = np.conj( self.real.data ) , np.conj( self.imag.data )
            indxp  = getMLvec( self.ntrunc , 'p' )
            indxpn = getMLvec( self.ntrunc , 'pn' )
            data = np.zeros( ( len(indxpn) , p.shape[1] ) , complex )
            for i,ml in enumerate(indxpn):
                m , l = ml[0] , ml[1]
                if m >= 0 :
                    k = indxp.index( ml )
                    data[i,:] = (-1)**m * sqrt2pi * ( p[k,:] + 1j*q[k,:] )
                else:
                    ml = (-m,l)
                    k = indxp.index( ml )
                    data[i,:] =  sqrt2pi * ( np.conj( p[k,:] ) + 1j*np.conj( q[k,:] ) )
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
    The class is designed to evaluate expression () in [paper/note name].  The various terms in this expression come from the different input objects in the constructor.  Their frequencies are turned into a common set, and the expression is evaluated at these frequencies.

    """
    def __init__(self, OrfMultipleMoments , PSDDIC , TimeSeries , 
                 GwSpectralSlope ,
                 day , cIJdir='./' , XXXdir='./' ):
        """
        Create an instance of the class XXX
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
#        self.sourcename = sourcename
        self.day = day
        self.cIJdir = cIJdir
        self.XXXdir = XXXdir
######################################################################
        return

    def getCrossSpectraData(self,window='None',fftlength=None):
        """Return cross spectral data

        Keyword parameters:
        window -- amplitude of window in numpy array
        fftlength -- length of fast fourier transform (fft)

        Output:
        f -- frequency
        C -- cross spectral data

        """
        if window == 'None' :
            window = np.ones( self.ts.t.data.shape[0] )
        elif window == 'hanning' :
            window = sp.hanning( len( self.ts.t.data ) )
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
        fcoarse = coarsefrequency( self.orf.f , self.psd['f'] , self.fcsdata )
        self.fcoarse = fcoarse
        return fcoarse

    def getintegrand(self,flow=1.4e-4,fhigh=1e-1 , lmax=20 , window='None' ):
        """Return the integrand in the expression for the dirty map

        Keyword parameter:
        flow -- lower frequency limit (default 1.4e-4Hz)
        fhigh -- higher frequency limit (default 1e-1Hz)
        lmax -- maximum l
        window -- window used on time-series before being ffted
        
        Output:
        f -- frequencies
        integrand -- the integrand

        Note:
        See expression () in [paper/notes name]
        
        """
        print 'Calculating XXX...'
        self.getCrossSpectraData( window = window , fftlength = None )
        self.getCoarseFreqs()

        if  fhigh < flow :
            print "fhigh must be >= flow"
            raise ValueError
        elif flow < self.fcoarse.data[0] :
            print "flow must be >= %f" % self.fcoarse.data[0]
            raise ValueError
        elif fhigh > self.fcoarse.data[-1] :
            print "fhigh must be <= %f" % self.fcoarse.data[-1]
            raise ValueError

        glm = self.orf.getMultipleMoments( 'pn' , lmax )
        pII = self.psd[ 2*self.pair[0] ]
        pJJ = self.psd[ 2*self.pair[1] ]
        csdata = self.csdata
        
###################### save csdata , pII , pJJ ########################(2)
#        cIJdict = {}
#        cIJdict['csdata']  = csdata
#        cIJdict['pII']     = pII
#        cIJdict['pJJ']     = pJJ
#        
#        file = open( self.cIJdir + self.sourcename + '_obs-s%03d-cIJ.pkl' % 
#                     self.day , 'wb' )
#        cpkl.dump( cIJdict , file , -1 )
#        file.close()
#########################################################################
        
        inputs = { 'Offset1'  : self.fcoarse.Offset1 ,
                   'Cadence1' : self.fcoarse.Cadence1 ,
                   'N1'       : len( self.fcoarse.data ) }
        glm = glm.coarsegrain( **inputs )
        pII , pJJ  = pII.coarsegrain( **inputs ) , pJJ.coarsegrain( **inputs )
        csdata = csdata.coarsegrain( **inputs )

######################## save csdata , pII , PJJ ####################(1)
        cIJdict = { 'f':self.fcoarse , 2*self.pair[0]:pII , 2*self.pair[1]:pJJ , self.pair:csdata }
        if self.cIJdir not in glob.glob( self.cIJdir ) :
            os.system( 'mkdir -p %s' % self.cIJdir )
        file = open( self.cIJdir + '/d%03d.pkl' % self.day , 'wb' )
        cpkl.dump( cIJdict , file , -1 ) ; file.close()
###########################################################################

        H = self.fcoarse.data**self.alpha

########## Estimator containing only +ve frequencies ######################
#        HoPP = H / ( np.real( pII.data ) * np.real( pJJ.data ) )
#        data = np.conj( glm.data ) * HoPP * csdata.data
###########################################################################
########## Estimator containing +ve and -ve frequencies ###################
        datap = np.conj( glm.data )
        HoPP = H / ( np.real( pII.data ) * np.real( pJJ.data ) )

        datan = np.zeros( datap.shape , dtype='complex' )
        indxpn = getMLvec( lmax , 'pn' )
        for i,ml in enumerate( indxpn ) :
            m , l = ml
            k = indxpn.index( ( -m , l ) )
            datan[ i ] = (-1)**m * glm.data[ k ]
                        
        data = ( datap * csdata.data + datan * np.conj( csdata.data ) ) * HoPP
###########################################################################
        
        ilow  = np.round( (flow - self.fcoarse.Offset1) / self.fcoarse.Cadence1 )
        ihigh = np.round( (fhigh - self.fcoarse.Offset1) / self.fcoarse.Cadence1 )
        data = data[ : , ilow : ihigh + 1 ]
        fdata = self.fcoarse.data[ ilow : ihigh + 1 ]
        fscale = {'Offset1' : fdata[0] , 
                  'Cadence1': fdata[1] - fdata[0] }
        f         = Coarsable( fdata , **fscale ) 
        integrand = Coarsable( data  , **fscale )
        return f , integrand
    def getsummand(self,flow=1.4e-4,fhigh=1e-1 , lmax=20 , window = 'None' ):
        """Return the summand in the expression for the dirty map

        Keyword parameters:
        flow -- lower frequency limit (default 1.4e-4Hz)
        fhigh -- higher frequency limit (default 1e-1Hz)
        lmax -- maximum l
        window -- window applied to the time-series before it is fft'ed.

        Output:
        summand -- the summand

        Note:
        See expression () in '[paper/note name]'.  The summand returned here is the daily dirty map.
        
        """
        print 'Calculating XX...'
        f , integrand = self.getintegrand( flow , fhigh , lmax=lmax , window = window )
        data = np.sum( integrand.data , 1 )
        summand = Coarsable( data )
        print 'done'
        return summand



class XXX_test(object):

    def __init__(self, OrfMultipleMoments , PSDDIC , CSDDIC , 
                 GwSpectralSlope ,
                 day , cIJdir='./' , XXXdir='./' ):
        self.csd     = CSDDIC
        self.psd     = PSDDIC
        self.orf     = OrfMultipleMoments
        self.alpha   = GwSpectralSlope
        self.pair    = self.orf.Antenna

        self.day = day
        self.cIJdir = cIJdir
        self.XXXdir = XXXdir
        return

    def getCrossSpectraData(self,fftlength=None):
        f = self.csd[ 'f' ]
        C = self.csd[ self.pair ]
        self.fcsdata , self.csdata = f , C
        return f , C

    def getCoarseFreqs(self):
        self.getCrossSpectraData()
        fcoarse = coarsefrequency( self.orf.f , 
                                   self.psd['f'] , 
                                   self.fcsdata )
        self.fcoarse = fcoarse
        return fcoarse

    def getintegrand(self,flow=1.4e-4,fhigh=1e-1 , lmax=20 ):
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

        glm = self.orf.getMultipleMoments( 'pn' , lmax )
        pII = self.psd[ 2*self.pair[0] ]
        pJJ = self.psd[ 2*self.pair[1] ]
        csdata = self.csdata
        
        inputs = { 'Offset1'  : self.fcoarse.Offset1 ,
                   'Cadence1' : self.fcoarse.Cadence1 ,
                   'N1'       : len( self.fcoarse.data ) }
        glm = glm.coarsegrain( **inputs )
        pII , pJJ  = pII.coarsegrain( **inputs ) , pJJ.coarsegrain( **inputs )
        csdata = csdata.coarsegrain( **inputs )

######################## save csdata , pII , PJJ ####################(1)
        cIJdict = { 'f':self.fcoarse , 2*self.pair[0]:pII , 2*self.pair[1]:pJJ , self.pair:csdata }
        if self.cIJdir not in glob.glob( self.cIJdir ) :
            os.system( 'mkdir -p %s' % self.cIJdir )
        file = open( self.cIJdir + '/d%03d.pkl' % self.day , 'wb' )
        cpkl.dump( cIJdict , file , -1 ) ; file.close()

###########################################################################

        H = self.fcoarse.data**self.alpha

        datap = np.conj( glm.data )
        HoPP = H / ( np.real( pII.data ) * np.real( pJJ.data ) )

        datan = np.zeros( datap.shape , dtype='complex' )
        indxpn = getMLvec( lmax , 'pn' )
        for i,ml in enumerate( indxpn ) :
            m , l = ml
            k = indxpn.index( ( -m , l ) )
            datan[ i ] = (-1)**m * glm.data[ k ]
                        
        data = ( datap * csdata.data + datan * np.conj( csdata.data ) ) * HoPP
        
        ilow  = np.round( (flow - self.fcoarse.Offset1) / 
                             self.fcoarse.Cadence1 )
        ihigh = np.round( (fhigh - self.fcoarse.Offset1) / 
                             self.fcoarse.Cadence1 )
        data = data[ : , ilow : ihigh + 1 ]
        fdata = self.fcoarse.data[ ilow : ihigh + 1 ]
        fscale = {'Offset1' : fdata[0] , 
                  'Cadence1': fdata[1] - fdata[0] }
        f         = Coarsable( fdata , **fscale ) 
        integrand = Coarsable( data  , **fscale )
        return f , integrand

    def getsummand(self,flow=1.4e-4,fhigh=1e-1 , lmax=20 ):
        print 'Calculating XX...'
        f , integrand = self.getintegrand( flow , fhigh , lmax=lmax )
        data = np.sum( integrand.data , 1 )
        summand = Coarsable( data )
        print 'done'
        return summand



def get_covariance_bias_matrix_for_the_day( orf , psd , csd , alpha , day , flow , fhigh , lmax ) :
    """
    returns the 'bias matrix' for the covariance of the estimators
    INPUT:
    orf --- OrfMultipleMoments object
    psd --- dictionary of psds
    csd --- dictionary of csds
    alpha --- GW spectral slope
    day --- day
    flow --- lower limit of integration over frequency
    fhigh --- higher limit of integration over frequency
    lmax --- maximum l to include
    OUTPUT:
    biasm --- 'bias matrix'
    """
    IJ , glm = orf.Antenna , orf.getMultipleMoments( 'pn' , lmax )
    cIJ , pII , pJJ = csd[ IJ ] , psd[ 2*IJ[0] ] , psd[ 2*IJ[1] ]
    fcoarse = coarsefrequency( orf.f , psd['f'] , csd['f'] )

    inputs = { 'Offset1':fcoarse.Offset1 , 'Cadence1':fcoarse.Cadence1 , 'N1':len(fcoarse.data) }
    glm = glm.coarsegrain( **inputs )
    cIJ , pII , pJJ = cIJ.coarsegrain( **inputs ) , pII.coarsegrain( **inputs ) , pJJ.coarsegrain( **inputs )

    H = fcoarse.data**alpha
    H2oP2P2 = H**2 / ( np.real( pII.data )**2 * np.real( pJJ.data )**2 )

    indxpn = getMLvec( lmax , 'pn' )
    Nml = len( indxpn )
    Npml = int( .5 * ( Nml - ( lmax + 1 ) ) )

    glmmdata = np.zeros( glm.data.shape , dtype='complex' )
    glmmdata[ : lmax+1 ] = glm.data[ : lmax+1 ]
    glmmdata[ lmax+1 : lmax+1 + Npml ] = glm.data[ lmax+1 + Npml : ]
    glmmdata[ lmax+1 + Npml : ] = glm.data[ lmax+1 : lmax+1 + Npml ]
    minus1tom = np.zeros( ( Nml, ) )
    for i,ml in enumerate( indxpn ) :
        minus1tom[ i ] = ( -1 )**ml[0]

    minus1tomTglmmdata = np.transpose( np.transpose( glmmdata ) * minus1tom )

    ilow  = int( np.round( (flow - fcoarse.Offset1) / fcoarse.Cadence1 ) )
    ihigh = int( np.round( (fhigh - fcoarse.Offset1) / fcoarse.Cadence1 ) )
    biasm = np.zeros( ( Nml , Nml ) , dtype='complex' )
    for ii in range( ilow , ihigh+1 ) :
        p = glm.data[ : , ii ]
        q = minus1tomTglmmdata[ : , ii ]
        A = cIJ.data[ ii ]**2 * np.outer( np.conj( p ) , np.conj( q ) )
        dbiasm = H2oP2P2[ii] * ( A + np.conj( np.transpose( A ) ) )
        biasm += dbiasm

    print 'PSI shape' , biasm.shape
    print 'PSI sample' , biasm[ :3 , :3 ]
    return biasm 




class GGG(object):
    """GGG is for getting the ingredients for the Fisher matrix

    Attributes:
    

    Methods:

    
    """
    def __init__(self, OrfMultipleMoments , PSDDIC , GwSpectralSlope , day , cIIdir='./' ):
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

        self.cIIdir = cIIdir
        return
    def getCoarseFreqs(self):
        """Return frequencies for coarsegraining 

        Output:
        fcoarse -- frequencies in a Coarsable           

        Note:
        This takes not input parameters. It also assigns the coarsegrain frequencies as a new attirbute of the instance.         

        """
        fcoarse = coarsefrequency( self.orf.f , self.psd['f'] )
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
        pII = self.psd[ 2*self.pair[0] ]
        pJJ = self.psd[ 2*self.pair[1] ]
        inputs = { 'Offset1'  : self.fcoarse.Offset1 ,
                   'Cadence1' : self.fcoarse.Cadence1 ,
                   'N1'       : len( self.fcoarse.data ) }
        glm = glm.coarsegrain( **inputs )
        pII , pJJ = pII.coarsegrain( **inputs ) , pJJ.coarsegrain( **inputs )
        H = self.fcoarse.data**self.alpha

        ml1m , ml2m = ( -ml1[0] , ml1[-1] ) , ( -ml2[0] , ml2[-1] )
        indxpn = getMLvec( self.orf.ntrunc , 'pn' )
        k1 , k2 = indxpn.index( ml1 ) , indxpn.index( ml2 )
        k3 , k4 = indxpn.index( ml1m ) , indxpn.index( ml2m )
        
        data = ( np.conj( glm.data[k1,:] ) * glm.data[k2,:] + (-1)**( ml1[0] + ml2[0] ) * np.conj( glm.data[k4,:] ) * glm.data[k3,:] ) * H**2 / np.real( pII.data ) / np.real( pJJ.data )

        ilow  = int( np.round( (flow - self.fcoarse.Offset1) /
                                  self.fcoarse.Cadence1 ) )
        ihigh = int( np.round( (fhigh - self.fcoarse.Offset1) /
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
        
        glm = self.orf.getMultipleMoments( 'pn' , lmax )
        pII = self.psd[ 2*self.pair[0] ]
        pJJ = self.psd[ 2*self.pair[1] ]

        inputs = { 'Offset1'  : self.fcoarse.Offset1 ,
                   'Cadence1' : self.fcoarse.Cadence1 ,
                   'N1'       : len( self.fcoarse.data  ) }
        glm = glm.coarsegrain( **inputs )
        pII , pJJ = pII.coarsegrain( **inputs ) , pJJ.coarsegrain( **inputs )
######################### save pII , PJJ ###################################
#        cIIdict = { 'f':self.fcoarse , 2*self.pair[0]:pII , 2*self.pair[1]:pJJ }
#        if self.cIIdir not in glob.glob( self.cIIdir ) :
#            os.system( 'mkdir -p %s' % self.cIIdir )
#        file = open( self.cIIdir + '/d%03d.pkl' % self.day , 'wb' )
#        cpkl.dump( cIIdict , file , -1 ) ; file.close()
############################################################################
        H = self.fcoarse.data**self.alpha

        ilow  = int( np.round( (flow - self.fcoarse.Offset1) /
                                  self.fcoarse.Cadence1 ) )
        ihigh = int( np.round( (fhigh - self.fcoarse.Offset1) /
                                  self.fcoarse.Cadence1 ) )


        H2oPP = H**2 / ( np.real( pII.data ) * np.real( pJJ.data ) )

        indxpn = getMLvec( lmax , 'pn' )
        Nml = len( indxpn )

        Npml = int( .5 * ( Nml - ( lmax + 1 ) ) )
        glmmdata = np.zeros( glm.data.shape , complex )
        glmmdata[ : lmax+1 ] = glm.data[ : lmax+1 ]
        glmmdata[ lmax+1 : lmax+1 + Npml ] = glm.data[ lmax+1 + Npml : ]
        glmmdata[ lmax+1 + Npml : ] = glm.data[ lmax+1 : lmax+1 + Npml ]

        minus1tom = np.zeros( ( Nml, ) )
        for i,ml in enumerate( indxpn ) :
            minus1tom[ i ] = ( -1 )**ml[0]

        minus1tomTglmmdata = np.transpose( np.transpose( glmmdata ) * minus1tom )

        data = np.zeros( ( Nml , Nml ) , dtype='complex' )
#        "~~ (l,m,l',m')s to check below"
#        mlmls = [ ( (0,0),(0,0) ) , ( (-2,2),(1,5) ), ( (4,6),(-3,4) ) , ( (2,14),(-6,12) ) ]
#        mlmlsm = [ ( ( -mlml[0][0] , mlml[0][1] ) , ( -mlml[1][0] , mlml[1][1] ) ) for mlml in mlmls ]
#        inds  = [ ( indxpn.index( mlml[0] ) , indxpn.index( mlml[1] ) ) for mlml in mlmls ]
#        indsm = [ ( indxpn.index( mlml[0] ) , indxpn.index( mlml[1] ) ) for mlml in mlmlsm ]
#        Nf = ihigh - ilow + 1
#        cheshs = [ np.zeros( (Nf,) , dtype='complex' ) for k in range( len( mlmls ) ) ]
#        "~~"
        for ii in range( ilow , ihigh+1 ) :
            p = np.conj( glm.data[ : , ii ] )
            q = minus1tomTglmmdata[ : , ii ]
            pp = np.outer( p , np.conj(p) )
            qq = np.outer( q , np.conj(q) )
            ddata = H2oPP[ii] * ( pp + qq )
            data += ddata
            
#            "~~ Get value for some (lm,l'm') at this frequency"
#            for k in range( len( mlmls ) ) :
#                cheshs[k][ ii - ilow ] = ddata[ inds[k] ]  
#            "~~"
            
#        "~~ Write H2oPP, cheshs and glms to disk"
#        cheshdirs = [ 'cheshs/chesh1/' , 'cheshs/chesh2/' , 'cheshs/chesh3/' , 'cheshs/chesh4/' ]
#        for cheshdir in cheshdirs :
#            if cheshdir not in glob.glob( cheshdir ) :
#                os.system( 'mkdir -p %s' % cheshdir )
#        
#        f_chesh = self.fcoarse.data[ ilow:ihigh+1 ]
#        cheshdicts = [ { 'f' : f_chesh ,
#                         'H2oPP' : H2oPP[ ilow:ihigh+1 ] ,
#                         'chesh' : cheshs[k] ,
#                         'ml1' : mlmls[k][0] , 'ml2' : mlmls[k][1] ,
#                         'glm1' : glm.data[ inds[k][0] , ilow:ihigh+1 ] ,
#                         'glm2' : glm.data[ inds[k][1] , ilow:ihigh+1 ] ,
#                         'glm1m' : glm.data[ indsm[k][0] , ilow:ihigh+1 ] ,
#                         'glm2m' : glm.data[ indsm[k][1] , ilow:ihigh+1 ]
#                         }
#                       for k in range( len(mlmls) ) ]
#        
#        for k in range( len( mlmls ) ) :
#            file = open( cheshdirs[k]+'chesh-s%03d.pkl' % self.day , 'wb' );
#            cpkl.dump( cheshdicts[k] , file , -1 )
#            file.close()
#        "~~"
        summand = Coarsable( data )

        print 'done'
        return summand





def Adjoint( matrix ) :
    return np.transpose( np.conj( matrix ) )


def IsPositiveDefinite( A , smallest=0 ) :
    """
    INPUT:
    A --- matrix ( numpy array )
    smallest --- print  the smallest 'smallest' Evalues
    OUTPUT:
    Boolean --- True or False
    """
    d , P = sp_linalg.eig( A )
    posdef = True
    for ii in range( d.shape[0] ) :
        if np.real( d[ ii ] ) <= 0 :
            posdef = False
            break
    if smallest :
        print np.sort( d )[ :smallest ]
    return posdef


def IsUeqV( A , abscomp=False ) :
    """
    INPUT:
    A --- matrix (numpy array)
    abscomp --- Boolean, if True returns (abs of mean) and (max abs of diff)
    OUTPUT:
    Boolean --- True or False
    """
    U , s , Vh = sp_linalg.svd( A )
    if abscomp :
        max_abs_diff = np.max( np.abs( U - Adjoint(Vh) ) )
        abs_mean = np.mean( np.abs( U ) )
        print "max. abs(diff):" , max_abs_diff
        print "mean of abs(matrix):" , abs_mean
    return np.allclose( U , Adjoint(Vh) )


def DiagonalsAllPositive( A , smallest=0 ) :
    """
    INPUT:
    A --- matrix (numpy array)
    OUTPUT:
    Boolean --- True or False
    """
    d = np.diag( A )
    allpos = True
    for ii in range( d.shape[0] ) :
        if np.real( d[ ii ] ) <= 0 :
            allpos = False
            break
    if smallest :
        print np.sort( d )[ :smallest ]
    return allpos


def ConditionNumber( A ) :
    """
    INPUT:
    A --- matrix ( numpy array )
    OUTPUT:
    cn --- condition number of A
    """
    U , s , Vh = sp_linalg.svd( A )
    cn = np.max(s) / np.min(s)
    return cn


def InvertWell( A , abscomp=False ) :
    """
    INPUT:
    A --- matrix ( numpy array )
    abscomp --- Boolean, if True returns (abs of mean) and (max abs of diff)
    OUTPUT:
    Boolean --- True or False
    """
    U , s , Vh = sp_linalg.svd( A )
    sinv = 1 / s
    Ainv = np.dot( Adjoint(Vh) , np.dot( np.diag(sinv) , Adjoint(U) ) )
    AinvA = np.dot( Ainv , A )
    I = np.identity( A.shape[0] )
    if abscomp :
        max_abs_diff = np.max( np.abs( AinvA - I ) )
        abs_mean = np.mean( np.abs( AinvA ) )
        print "max. abs(diff):" , max_abs_diff
        print "mean of abs(matrix):" , abs_mean
    return np.allclose( AinvA , I )

    
def IsHermitian( matrix , abscomp = False ) :
    """
    INPUT:
    matrix --- numpy array
    abscomp --- Boolean, if True returns (abs of mean) and (max abs of diff)
    OUTPUT:
    Boolean --- True or False
    """
    if abscomp :
        max_abs_diff = np.max( np.abs( matrix - Adjoint(matrix) ) )
        abs_mean = np.mean( np.abs( matrix ) )
        print "max. abs(diff):" , max_abs_diff
        print "mean of abs(matrix):" , abs_mean
    return np.allclose( Adjoint(matrix) , matrix )


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
                                                                                                        


def avg_psdPKL_avg(day,lpsdpath,rpsdpath,psddir):
    file = open( lpsdpath , 'rb' )
    lpsd = cpkl.load( file )
    file.close()
    file = open( rpsdpath , 'rb' )
    rpsd = cpkl.load( file )
    file.close()

    PAA = ( lpsd['AA'].data + rpsd['AA'].data ) / 2
    PEE = ( lpsd['EE'].data + rpsd['EE'].data ) / 2
    PTT = ( lpsd['TT'].data + rpsd['TT'].data ) / 2

    fdata = lpsd['f'].data ; f0 = fdata[0] ; df = fdata[1] - fdata[0]
    f   = Coarsable( fdata , Offset1=f0 , Cadence1=df ) 
    PAA = Coarsable( PAA , Offset1=f0 , Cadence1=df )
    PEE = Coarsable( PEE , Offset1=f0 , Cadence1=df )
    PTT = Coarsable( PTT , Offset1=f0 , Cadence1=df )

    psd = { 'f':f , 'AA':PAA , 'EE':PEE , 'TT':PTT }

#    lpsdname = os.path.basename( lpsdpath )
#    psdname = re.split( '-' , lpsdname , 1 )[0] + '-s%03d-psd.pkl' % day
    psdpath = psddir + 'd%03d.pkl' % day
    file = open( psdpath , 'wb' ) ; cpkl.dump( psd , file , -1 ) ; file.close()
    return



def IsMinusOneP( matrix , lmax ) :
    indxpn = getMLvec( lmax , 'pn' )
    Nml = len( indxpn )
    Npml = int( 0.5 * ( Nml - ( lmax + 1 ) ) )

    minus1tomPm = np.zeros( ( Nml , Nml ) )
    for i1,ml1 in enumerate( indxpn ) :
        for i2 , ml2 in enumerate( indxpn ) :
            minus1tomPm[ i1 , i2 ] = ( -1 )**( ml1[0] + ml2[0] )
    datamm = np.zeros( matrix.shape , dtype = matrix.dtype )
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
    return np.allclose( minus1tomPmTdatamm , np.conj( matrix ) )

    



class FisherMatrix(object):
    def __init__( self , fishpath , lmax = None ):
        file = open( fishpath , 'rb' ) ; fishdict = cpkl.load( file ) ; file.close()
        if lmax == None :
            self.ntrunc = fishdict['ntrunc']
            self.fish = fishdict['G'].data
        else :
            self.ntrunc = lmax
            self.fish = get_lmax_subset_from( fishdict['G'].data , lmax )

        if not IsHermitian( self.fish ) :
            print "Warning: Fisher matrix NOT Hermitian."

        if not IsMinusOneP( self.fish , self.ntrunc ) :
            print "Warning: Fisher matrix does not satisfy the minus 1 to the m+m' identity."

        self.decomposed  = False
        self.regularised = False
        self.N_keptSV = ( self.ntrunc + 1 )**2
        return


    def svd(self):
        self.U, self.s, self.Vh = sp_linalg.svd(self.fish)
        self.decomposed = True
        print 'Condition number = %f' % ( np.max( self.s ) / np.min( self.s ) )
        return self.U , self.s , self.Vh


    def regularise(self,regMethod=2,regCutoff=2./3):
        if not self.decomposed :
            self.svd()
        if regMethod in [ 1 , 11 ]:
            ind = np.where( self.s/self.s[0] >= regCutoff )[0]
            if len(ind) > 0 :
                ii = ind[-1]
            else:
                ii = 0
        elif regMethod in [ 2 , 22 ]:
            ii = int( np.round( len(self.s)*regCutoff ) ) - 1
        elif regMethod in [ 3 , 33 ]:
            ii = int( np.round( regCutoff ) ) - 1
        else:
            print 'Unknown regularisation method.'
            raise Exception
        if ii < 1 :
            ii = 0
            if ii > len(self.s)-1 :
                ii = len(self.s)-1 #-1
        self.regs = np.zeros(self.s.shape,self.s.dtype)
        if regMethod in [ 1 , 2 , 3 ]:
            self.regs[:ii+1] , self.regs[ii+1:] = self.s[:ii+1] , self.s[ii]
        elif regMethod in [ 11 , 22 , 33 ]:
            self.regs[:ii+1] , self.regs[ii+1:] = self.s[:ii+1] , np.inf
        else:
            print 'Unknown regularisation method.'
        self.regularised = True
        self.N_keptSV = ii + 1
        return self.regs


    def invert(self):
        if not self.decomposed:
            self.svd()
        self.sinv = 1. / self.s
        Sinv = np.diag( self.sinv )
        Uh = np.transpose(np.conj(self.U))
        V  = np.transpose(np.conj(self.Vh))
        fishinv = np.dot( V , np.dot( Sinv , Uh ) )
        for ii in range( fishinv.shape[0] ) :
            fishinv[ ii , ii ] = np.real( fishinv[ ii , ii ] )
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
        Sinv = np.diag( self.regsinv )
        Uh = np.transpose(np.conj(self.U))
        V  = np.transpose(np.conj(self.Vh))
        regfishinv = np.dot( V , np.dot( Sinv , Uh ) )

        for ii in range( regfishinv.shape[0] ) :
            regfishinv[ ii , ii ] = np.real( regfishinv[ ii , ii ] )
        self.regfishinv = regfishinv

        if not IsHermitian( self.regfishinv ) :
            print "Warning: Fisher matrix's inverse NOT Hermitian."
        if not IsMinusOneP( self.regfishinv , self.ntrunc ) :
            print "Warning: Fisher matrix's inverse does not satisfy the minus 1 to the m+m' identity."
        return self.regfishinv


    def getpCovar(self):
        if not self.regularised :
            self.regularise()
        Cov = np.diag( self.s / self.regs**2 )
        Uh = np.transpose(np.conj(self.U))
        V  = np.transpose(np.conj(self.Vh))
        self.pCovar = np.dot(V,np.dot(Cov,Uh))
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
    cleanmap = np.dot( fishinv , skymap.xlm )
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
    
    mlvec0 = getMLvec( lmax0 )
    mlvec = getMLvec( lmax )
    
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

    mlvec0 = getMLvec( lmax0 )

    indices = []
    for i,ml in enumerate( mlvec0 ) :
        if ml[-1] == l :
            indices += [i]
        else :
            continue
    return indices 
            
    
def get_lmax_subset_from( matrix , lmax , axis='All' ) :
    """
    returns a subset of matrix restricted to lmax
    INPUT:
    matrix --- 1/2-D numpy array whose side or sides correspond to (m,l)s
    lmax --- maximum l to restrict matrix to
    axis --- along which access to restric to lmax
    OUTPUT :
    new_A --- matrix, restricted to lmax
    """
    if axis == 0 :
        A = np.copy( matrix )
        if np.mod( np.sqrt( A.shape[0] ) - 1 , 1 ) != 0 :
            raise Exception , "Side must have length corresponding to an integer lmax."
        else:
            lmax0 = int( np.sqrt( A.shape[0] ) - 1 )
        
        indices = get_MLvec_indices( lmax0 , lmax )
        new_A = A[ indices ]
        return new_A

    elif axis == 1 :
        A = np.transpose( np.copy( matrix ) )
        new_A = np.transpose( get_lmax_subset_from( A , lmax , axis=0 ) )
        return new_A

    elif axis == 'All' :
        A = get_lmax_subset_from( np.copy( matrix ) , lmax , axis=0 )
        new_A = get_lmax_subset_from( A , lmax , axis=1 )
        return new_A
    
    else:
        raise ValueError , "axis must either be 0, 1, or 'All'"


    
def get_l_subset_from( matrix , l , axis='All' ) :
    """
    returns a subset of matrix restricted to l
    INPUT:
    matrix --- 1/2-D numpy array whose side or sides correspond to (m,l)s
    l --- l to restrict matrix to
    axis --- along which access to restric to l
    OUTPUT :
    new_A --- matrix, restricted to l
    """
    
    if axis == 0 :
        A = np.copy( matrix )
        if np.mod( np.sqrt( A.shape[0] ) - 1 , 1 ) != 0 :
            raise Exception , "Side must have length corresponding to an integer lmax."
        else:
            lmax0 = int( np.sqrt( A.shape[0] ) - 1 )

        indices = get_MLvec_indices_from_l( lmax0 , l )
        new_A = A[ indices ]
        return new_A

    elif axis == 1 :
        A = np.transpose( np.copy( matrix ) )
        new_A = np.transpose( get_l_subset_from( A , l , axis=0 ) )
        return new_A

    elif axis == 'All' :
        A = get_l_subset_from( np.copy( matrix ) , l , axis=0 )
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
    lmax = int( np.sqrt( covar.shape[0] ) - 1 )
    U , lats , lons = checkPIXELCONVERSION( lmax , nlat , nlon )
    Uh = np.transpose( np.conj( U ) )
    covarP = np.dot( U , np.dot( covar , Uh ) )
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
    sigmaP = np.sqrt( np.real( np.diag( covarP ) ) )
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
    lmax = np.sqrt( Nml ) - 1
    if ( lmax % 1 ) != 0 :
        raise Exception , 'Number of rows does not correspond to an integer lmax.'
    else:
        lmax = int( lmax )
        
    indxpn = getMLvec( lmax )

    ml_matrixl = list( ml_matrix )
    lm_matrixl = []
    for l in range( lmax+1 ) :
        for m in range( -l , l+1 ) :
            lm_matrixl += [ ml_matrixl[ indxpn.index( (m,l) ) ] ]
    lm_matrix = np.array( lm_matrixl )
    return lm_matrix


def get_Angular_PSD( mla ) :
    """
    Returns angular PSD from an array whose elements are arranged according to (m,l) from getMLvec()
    INPUT:
    mla --- multipole moments aranged by LISAresonse.getMLvec( msign='pn' ). (Contains m<0)
    OUTPUT:
    C --- angular power spectral density
    """
    Nml = mla.shape[0]
    lmax = np.sqrt( Nml ) - 1
    if ( lmax % 1 ) != 0 :
        raise Exception , 'Number of rows does not correspond to an integer lmax.'
    else:
        lmax = int( lmax )

    indxpn = getMLvec( lmax )

    mlal = list( mla )

    C = np.zeros( ( lmax+1 , ) )
    for l in range( lmax+1 ) :
        for m in range( -l , l+1 ) :
            mlv = mlal[ indxpn.index( (m,l) ) ]
            C[l] += np.abs( mlv )**2
        C[l] = C[l] / ( 2*l + 1 ) 
    return C

#Utilities4.py body
class SkyMap(object):

    def __init__( self ):
        self.__ntrunc = []
        self.__plm    = []
        self.__qlm    = []
        self.__xlm    = []

        self.__sky    = []
        self.__P      = []
        self.__Q      = []
        self.__X      = []

    def get_ntrunc( self ):
        return self.__ntrunc
    ntrunc = property( get_ntrunc )
    
    def get_plm( self ):
        return self.__plm
    plm = property( get_plm )
    
    def get_qlm( self ):
        return self.__qlm
    qlm = property( get_qlm )
    
    def get_xlm( self ):
        return self.__xlm
    xlm = property( get_xlm )

    def get_sky( self ):
        return self.__sky
    sky = property( get_sky )
    
    def get_P( self ):
        return self.__P
    P = property( get_P )

    def get_Q( self ):
        return self.__Q
    Q = property( get_Q )

    def get_X( self ):
        return self.__X
    X = property( get_X )

    def TotalPower( self , realimag='real' ):
        if realimag == 'real' :
            powerdensity = np.copy( self.P )
        elif realimag == 'imag' :
            powerdensity = np.copy( self.Q )
        else:
            raise Exception , "realimag must either be 'real' or 'imag'. "

        dthe = np.radians( self.sky.lats[0] - self.sky.lats[1] )
        dphi = np.radians( self.sky.lons[1] - self.sky.lons[0] )

        power = 0
        for i,elat in enumerate( list( self.sky.lats ) ) :
            for j,elon in enumerate( list( self.sky.lons ) ) :
                the = np.radians( 90 - elat )
                power += np.sin( the ) * powerdensity[ i,j ]**2
        power = dthe * dphi * power 
        return power


    def AngularPower( self , l=0 , realimag='real' ) :
        if realimag == 'real' :
            xlm = np.copy( self.plm )
        elif realimag == 'imag' :
            xlm = np.copy( self.qlm )
        else :
            raise Exception , "realimag must either be 'real' or 'imag'."

        mlvec = getMLvec( self.ntrunc )        

        angpower = 0
        for m in range( l+1 ) :
            dangpower = abs( xlm[ mlvec.index( (m,l) ) ] )**2
            if m == 0 :
                angpower += dangpower
            else :
                angpower += 2 * dangpower
        return angpower
            
    def FractionalPower( self , l=0 , realimag='real' ) :
        return self.AngularPower( l=l , realimag=realimag ) / self.TotalPower( realimag=realimag )



class xlmSkyMap( SkyMap ):

    def __init__( self , ntrunc=5 , xlm=None ):
        super( xlmSkyMap , self ).__init__()
        
        if xlm == None:
            self._SkyMap__xlm = np.zeros( ( (ntrunc+1)**2, ) , dtype='complex' )
            self._SkyMap__ntrunc = ntrunc
        else:
            ntrunc = np.sqrt( xlm.shape[0] ) - 1
            if ( ntrunc % 1 ) != 0 :
                print "The length of the entered numpy array does not correspond to an integer ntrunc."
                raise ValueError
            else:
                self._SkyMap__xlm =  np.copy( xlm )
                self._SkyMap__ntrunc = int( ntrunc )
        return

    def alter_xlm( self , add_on=False , *xlms ):
        if xlms == ():
            print "No arrays to add/overwrite."
            return
        
        for xlm in xlms:
            if xlm.shape != self.xlm.shape:
                continue
            if add_on == False:
                self._SkyMap__xlm =  np.copy( xlm )
            elif add_on == True:
                self._SkyMap__xlm += xlm 
        return

    def alter_ml( self , add_on=False , *mlvs ):
        if mlvs == ():
            print "No (m,l) to add/overwrite."
            return

        mlvec = getMLvec( self.ntrunc )

        for mlv in mlvs:
            m , l , value = mlv

            if (m,l) not in mlvec:
                print "(%d,%d) is out of range. lmax = %d . " % ( m , l ,self.ntrunc )
                continue

            indx = mlvec.index( (m,l) )
            if add_on==False:
                self._SkyMap__xlm[ indx ] = value
            elif add_on==True:
                self._SkyMap__xlm[ indx ] = self._SkyMap__xlm[ indx ] + value
        return

    def xlm_to_plm( self ):
        indxpn = getMLvec( self.ntrunc , 'pn' )
        plm = np.zeros( self.xlm.shape , dtype='complex' )
        for i,ml in enumerate( indxpn ) :
            m , l = ml[0] , ml[1]
            k = indxpn.index( (-m,l) )
            plm[i] = ( self.xlm[i] + (-1)**m * np.conj( self.xlm[k] ) ) / 2
        indxp = getMLvec( self.ntrunc , 'p' )
        self._SkyMap__plm = np.copy( plm[ : len(indxp) ] )
        return

    def xlm_to_qlm( self ):
        indxpn = getMLvec( self.ntrunc , 'pn' )
        qlm = np.zeros( self.xlm.shape , dtype='complex' )
        for i,ml in enumerate( indxpn ) :
            m , l = ml[0] , ml[1]
            k = indxpn.index( (-m,l) )
            qlm[i] = ( self.xlm[i] - (-1)**m * np.conj( self.xlm[k] ) ) / 2j
        indxp = getMLvec( self.ntrunc , 'p' )
        self._SkyMap__qlm = np.copy( qlm[ :len(indxp) ] )
        return

    def create_sky( self , nlon=6 , nlat=4 ):
        self._SkyMap__sky = mySpharmt( nlon , nlat )
        return
    
    def plm_to_P( self ):
        indxp = getMLvec( self.ntrunc , 'p' )
        to_pyspharm = np.zeros( ( len(indxp), ) )
        for i,ml in enumerate( indxp ):
            to_pyspharm[i] = (-1)**ml[0] / np.sqrt( 2*np.pi )

        self._SkyMap__P = self.sky.spectogrd( to_pyspharm * self.plm )
        return

    def qlm_to_Q( self ):
        indxp = getMLvec( self.ntrunc , 'p' )
        to_pyspharm = np.zeros( ( len(indxp), ) )
        for i,ml in enumerate( indxp ):
            to_pyspharm[i] = (-1)**ml[0] / np.sqrt( 2*np.pi )

        self._SkyMap__Q = self.sky.spectogrd( to_pyspharm * self.qlm )
        return

    def PQ_to_X( self ):
        if [] in [ self.P , self.Q ]:
            print "Make sure P and Q are calculated first."
            raise
        self._SkyMap__X = self.P + 1j*self.Q
        return


class XSkyMap( SkyMap ):
    
    def __init__( self , nlat=3 , nlon=4  , X=None ):
        super( XSkyMap , self ).__init__()

        if X == None:
            if nlat < 3 or nlon < 4:
                print "Make sure that nlat >=3, and nlon >=4."
                raise ValueError
            else:
                self._SkyMap__X   = np.zeros( ( nlat , nlon ) , dtype='complex' )
                self._SkyMap__sky = mySpharmt( nlon , nlat )
        else:
            nlat , nlon = X.shape
            self._SkyMap__X   = np.copy( X )
            self._SkyMap__sky = mySpharmt( nlon , nlat )
        return

    def alter_X( self , add_on=False , *Xs ):
        if Xs == ():
            print "No arrays to add/overwrite."
            return

        for X in Xs:
            if X.shape != self.X.shape:
                continue
            if add_on == False:
                self._SkyMap__X = np.copy( X )
            elif add_on == True:
                self._SkyMap__X += X
        return


    def X_to_P( self ):
        self._SkyMap__P = np.real( self.X )
        return

    def X_to_Q( self ):
        self._SkyMap__Q = np.imag( self.X )
        return

    def ntrunc_equals( self , ntrunc ):
        if ntrunc > self.sky.nlat - 1:
            print "Make sure ntrunc <= nlat - 1."
            raise ValueError
        else:
            self._SkyMap__ntrunc = ntrunc
        return

    def P_to_plm( self ):
        indxp = getMLvec( self.ntrunc , 'p' )
        to_standard = np.zeros( ( len(indxp), ) )

        for i, ml in enumerate( indxp ):
            to_standard[i] = (-1)**ml[0] * np.sqrt( 2*np.pi )

        self._SkyMap__plm = to_standard * self.sky.grdtospec( self.P , self.ntrunc )
        return

    def Q_to_qlm( self ):
        indxp =getMLvec( self.ntrunc , 'p' )
        to_standard = np.zeros( ( len(indxp), ) )

        for i, ml in enumerate( indxp ):
            to_standard[i] = (-1)**ml[0] * np.sqrt( 2*np.pi )
                        
        self._SkyMap__qlm = to_standard * self.sky.grdtospec( self.Q , self.ntrunc )
        return

    def plmqlm_to_xlm( self ):
        indxp  = getMLvec( self.ntrunc , 'p' )
        indxpn = getMLvec( self.ntrunc , 'pn' )
        xlm = np.zeros( ( len(indxpn), ) , dtype='complex' )
        for i,ml in enumerate( indxpn ):
            m , l = ml
            k = indxp.index( ( abs(m),l ) )
            if m >= 0:
                xlm[i] = self.plm[k] + 1j*self.qlm[k]
            else:
                xlm[i] = (-1)**m * ( np.conj( self.plm[k] ) + 1j*np.conj( self.qlm[k] ) )
        self._SkyMap__xlm = np.copy( xlm )
        return
        


#def Convolve( orf , skymap , GWSpectralSlope=-3 ):
#
#    if not orf.ntrunc == skymap.ntrunc:
#        print "Make sure that lmax of the orf and the skymap are equal."
#        raise ValueError
#    else:
#        g = orf.getMultipleMoments( 'pn' )
#        f = g.Offset1 + g.Cadence1 * np.arange( g.data.shape[1] ) 
#        H = f**GWSpectralSlope
#        data = H * np.sum( np.transpose( g.data ) * skymap.xlm , 1 )
#        cspectrum = Coarsable( data ,
#                                          Offset1 = g.Offset1, Cadence1 = g.Cadence1 )
#
#    return cspectrum
    
def Convolve( orf , skymap , GWSpectralSlope=-3 ) :
    lmax = min( orf.ntrunc , skymap.ntrunc )
    g = orf.getMultipleMoments( 'pn' , lmax )
    f = g.Offset1 + g.Cadence1 * np.arange( g.data.shape[1] )
    H = f**GWSpectralSlope
    glm = g.data
    xlm = get_lmax_subset_from( skymap.xlm , lmax , axis=0 )
    data = H * np.sum( np.transpose( glm ) * xlm , 1 )
    cspectrum = Coarsable( data , Offset1 = g.Offset1 , Cadence1 = g.Cadence1 )
    return cspectrum


class TheZs( object ):
    def __init__( self , N=10 , seed=None ):

        np.random.seed( seed )
        p = [ np.random.standard_normal( (N,) ) for i in range(3) ]
        q = [ np.random.standard_normal( (N,) ) for i in range(3) ]

        self.__z1 = ( p[0] + 1j*q[0] ) / np.sqrt( 2 )
        self.__z2 = ( p[1] + 1j*q[1] ) / np.sqrt( 2 )
        self.__z3 = ( p[2] + 1j*q[2] ) / np.sqrt( 2 )

        self.__seed = seed
        return

    def get_z1( self ):
        return self.__z1
    z1 = property( get_z1 )

    def get_z2( self ):
        return self.__z2
    z2 = property( get_z2 )

    def get_z3( self ):
        return self.__z3
    z3 = property( get_z3 )

    def get_seed( self ):
        return self.__seed
    seed = property( get_seed )



class ShortTermFT( object ):
    def __init__( self , s1=[] , s2=[] , s3=[] ):

        self.__s1 , self.__s2 , self.__s3 = s1 , s2 , s3

        return

    def get_s1( self ):
        return self.__s1
    s1 = property( get_s1 )

    def get_s2( self ):
        return self.__s2
    s2 = property( get_s2 )

    def get_s3( self ):
        return self.__s3
    s3 = property( get_s3 )

    def alter_s1( self , s1=[] ):
        self.__s1 = s1
        return

    def alter_s2( self , s2=[] ):
        self.__s2 = s2
        return

    def alter_s3( self , s3=[] ):
        self.__s3 = s3
        return



def CSpectra_to_ShortTermFT( cspecdict , seed=None ):

    cspecnames = [ '11' , '12' , '13' ,
                   '22' , '23' ,
                   '33' ]

    cspeclist = [ cspecdict[ cspecname ] for cspecname in cspecnames ]

    freqlist = [ ( cspec.Offset1 , cspec.Cadence1 , cspec.data.shape[0] ) for cspec in cspeclist ]
    for i in range( 1 , len( freqlist ) ):
        if not freqlist[0] == freqlist[i]:
            print "Make sure all cross-spectra's frequencies are the same."
            raise ValueError
    
    Nf       = cspeclist[0].data.shape[0]
    fOffset  = cspeclist[0].Offset1
    fCadence = cspeclist[0].Cadence1

    Z = TheZs( Nf , seed )

    sdatas = [ np.zeros( (Nf,) , dtype='complex' ) for i in range(3) ]

    for k in range( Nf ):

        clist = [ cspec.data[k] for cspec in cspeclist ]
        cmatrix = np.array( [ [ clist[0]               , clist[1]               , clist[2] ] ,
                                 [ np.conj( clist[1] ) , clist[3] , clist[4] ]  ,
                                 [ np.conj( clist[2] ) , np.conj( clist[4] ) , clist[5] ]
                                 ] )

        w , V = sp_linalg.eig( cmatrix )

        for ww in list( w ):
            if ww <= 0:
                print "Warning: trying to simulate an unphysical signal."

        Vadjoint = np.transpose( np.conjugate( V ) )
        if not np.allclose( np.dot( V , Vadjoint ) , np.diag( np.ones( V.shape[0] ) ) ):
            print "Warning: the matrix of eigenvectors is not unitary."
        
        s = ( Z.z1[k]*cmath.sqrt( np.real( w[0] ) )*V[:,0] + 
              Z.z2[k]*cmath.sqrt( np.real( w[1] ) )*V[:,1] +
              Z.z3[k]*cmath.sqrt( np.real( w[2] ) )*V[:,2]
              )


        sdatas[0][k] , sdatas[1][k] , sdatas[2][k] = np.conj( s )


    scoarsables = [ Coarsable( sdata , Offset1=fOffset , Cadence1=fCadence ) for sdata in sdatas ]    
    shortft = ShortTermFT( s1=scoarsables[0] , s2=scoarsables[1] , s3=scoarsables[2] )

    return shortft


class sTimeSeries( object ):
    def __init__( self , s1=[] , s2=[] , s3=[] ):

        self.__s1 , self.__s2 , self.__s3 = s1 , s2 , s3

        return

    def get_s1( self ):
        return self.__s1
    s1 = property( get_s1 )
    
    def get_s2( self ):
        return self.__s2
    s2 = property( get_s2 )
    
    def get_s3( self ):
        return self.__s3
    s3 = property( get_s3 )
    
    def alter_s1( self , s1=[] ):
        self.__s1 = s1
        return
    
    def alter_s2( self , s2=[] ):
        self.__s2 = s2
        return
    
    def alter_s3( self , s3=[] ):
        self.__s3 = s3
        return
    



def InverseFT( shortft , Neven=False , fNyq=False , tOffset=0 ):
    
    ftlist = [ shortft.s1 , shortft.s2 , shortft.s3 ]

    tslist = []
    for ft in ftlist:

        if ft.Offset1 == 0:
            pftdata = ft.data
        else:
            pftdata = np.array( [0] + list( ft.data ) )


        if Neven == False:
            nftdata = np.conj( np.flipud( pftdata ) )[:-1]
            ftdata = np.concatenate( ( pftdata , nftdata ) )

        elif ( Neven , fNyq ) == ( True , True  ):
            nftdata = np.conj( np.flipud( pftdata ) )[1:-1]
            ftdata = np.concatenate( ( pftdata , nftdata ) )

        elif ( Neven , fNyq ) == ( True , False ):
            nftdata = np.conj( np.flipud( pftdata ) )[:-1]
            ftdata = np.concatenate( ( pftdata , np.array( [0] ) , nftdata ) )

        N = ftdata.shape[0] 
        dt = 1. / ( N * ft.Cadence1 )
        norm = np.sqrt( N / ( 2.*dt ) )

        tsdata = sp.ifft( norm * ftdata  )
        tsdata = np.real( tsdata )
        tslist += [ Coarsable( data=tsdata , Offset1=tOffset , Cadence1=dt ) ]

    tsdict = {}
    tsdict['s1'] , tsdict['s2'] , tsdict['s3'] = tslist[0] , tslist[1] , tslist[2]

    return sTimeSeries( **tsdict )


def WindowAndJoin( lts , rts , leftover=None ):

    if lts.data.shape != rts.data.shape :
        raise Exception , "Left and right arrays' lengths not equal."
    else:
        N = lts.data.shape[0]

    if lts.Cadence1 != rts.Cadence1 :
        raise Exception , "Left and right cadences not equal."
  
    if N % 2 == 0:
        lmid = N / 2
    else:
        lmid = ( N - 1 ) / 2

    if not ( lts.Cadence1 * lmid  + lts.Offset1 ) == rts.Offset1 :
        raise Exception , "Left's midpoint's offset not equal to right's offset."

    if not leftover == None :
        if not leftover.Cadence1 == lts.Cadence1:
            raise Exception , "Input leftover's cadence not equal to left's cadence."
        if not leftover.Offset1 == lts.Offset1 :
            raise Exception , "Input leftover's offset not equal to left's offset."
        if not leftover.data.shape  == lts.data[ :lmid ].shape :
            raise Exception , "Input leftover's length not equal to left's left half's length."

        window = np.sin( np.pi/N * np.arange( N ) )
        winlts = window * lts.data
        winrts = window * rts.data
        
        jointdata = np.concatenate( ( leftover.data + winlts[ :lmid ] , winlts[ lmid: ] + winrts[ :N-lmid ] ) )
        joint = Coarsable( jointdata , Offset1=lts.Offset1 , Cadence1=lts.Cadence1 )
        
    else:

        window = np.sin( np.pi/N * np.arange( N ) )
        winrts = window * rts.data

        jointdata = np.concatenate( ( lts.data[ :lmid ] ,
                                         window[ lmid: ] * lts.data[ lmid: ] + winrts[ :N-lmid ] )
                                       )
        joint = Coarsable( jointdata , Offset1=lts.Offset1 , Cadence1=lts.Cadence1 )
        
    whatsleftdata = winrts[ N-lmid: ]
    whatsleft = Coarsable( whatsleftdata ,
                                      Offset1 = rts.Offset1 + ( N - lmid) * lts.Cadence1 ,
                                      Cadence1 = lts.Cadence1 )
    
    return joint , whatsleft




def checkPIXELCONVERSION( lmax , nlat , nlon ) :
    """
    INPUT:
    lmax -- maximum degree of SpH, l
    nlat -- number of latitudes
    nlon -- number of longitudes 
    OUTPUT:
    U -- matrix for converting a matrix from SpH to PIXELS
    lats -- latitudes, repeated nlon times
    lons -- longitudes, repeated nlat times
    """
    indxpn = getMLvec( lmax , 'pn' )
    U = np.empty( ( nlat*nlon , len(indxpn) ) , dtype='complex' )
    for i,ml in enumerate( indxpn ) :
        M = xlmSkyMap( ntrunc=ml[-1] )
        M.alter_ml( False , ( ml[0] , ml[-1] , 1 ) )
        M.xlm_to_plm()
        M.xlm_to_qlm()
        M.create_sky( nlat=nlat , nlon=nlon )
        M.plm_to_P()
        M.qlm_to_Q()
        M.PQ_to_X()
        U[ :,i ] = np.reshape( M.X , ( nlat*nlon, ) )

    lats = np.outer( M.sky.lats , np.ones( nlon ) )
    lons = np.outer( np.ones( nlat ) , M.sky.lons )
    lats = np.reshape( lats , ( nlat*nlon , ) )
    lons = np.reshape( lons , ( nlat*nlon , ) )
    return U , lats , lons


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
    lmax = int( np.sqrt( covar.shape[0] ) - 1 )
    U , lats , lons = checkPIXELCONVERSION( lmax , nlat , nlon )
    Uh = np.transpose( np.conj( U ) )
    covarP = np.dot( U , np.dot( covar , Uh ) )
    return covarP , lats , lons
