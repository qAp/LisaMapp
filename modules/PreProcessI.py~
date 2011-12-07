import lisaxml
import AnisotropySearch as AS
import matplotlib; matplotlib.use( 'Agg' )
import pylab as pl




class TimeSeries(object):
    """TimeSeries class handles time-series in the analysis

    Attributes:
    t -- time 
    A -- time-series A
    E -- time-series E
    T -- time-series T (all the above 4 attributes are AS.Coarsables)

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
        The PSD is calculated using pl.psd() with an additional factor of dt in the normalisation.
        This gives the correct dimension.

        """
        x = getattr(self,xstr)
        fs = 1. / x.Cadence1
        P , f = pl.psd( x.data , nfft , fs , noverlap=noverlap , scale_by_freq=False )
        P = P / fs
        if coarsable:
            scale = {'Offset1' : f[0] ,
                     'Cadence1': f[1] - f[0] }
            f = AS.Coarsable( f , **scale )
            P = AS.Coarsable( P , **scale )
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
        The CPSD is calculated using pl.csd() with an additional factor of dt in the normalisation.
        If xstr='A' and ystr='E', then the CSD between A and E is returned.

        """
        x = getattr(self,xstr)
        y = getattr(self,ystr)
        fs = 1. / x.Cadence1
        P , f = pl.csd( x.data , y.data , nfft , fs , noverlap=noverlap )
        P = P / fs
        if coarsable:
            scale = {'Offset1' : f[0] ,
                     'Cadence1': f[1] - f[0] }
            f = AS.Coarsable( f , **scale )
            P = AS.Coarsable( P , **scale )
        return f , P


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
    cpkl.dump( psddict , file , -1 )
    file.close()
    return



