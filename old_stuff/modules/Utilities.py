from pylab import *
import synthlisa
import lisaxml
import numpy
import time
import re
import glob
import os
import sys


class TimeSeries(object):
    def __init__(self):
        self.data = []
        self.deltat = []
        return
    def setfbase(self,fbase):
        self.fbase = fbase
        return
    def setphase(self,phase):
        self.phase = phase
        return


class FreqSeries(object):
    def __init__(self):
        self.data = []
        self.flow = []
        self.deltaf = []
        self.symmetry = []
        return


def sumTerms(x,jlow,jhigh):
    Ny = len(jlow)
    y = zeros(Ny,complex64)
    for k in range(Ny):
        y[k] = sum(x[jlow[k]:jhigh[k]])
    return y


def coarsegrain(x,flowy,deltafy,Ny):
    """
    returns a frequency-series of a coarser resolution from the given one.
    INPUT:
    x - a FreqSeries object
    flowy - lower frequency of the coarser series
    deltafy - frequency step of the coarser series
    Ny - number of discrete frequencies of the coarser series
    OUTPUT:
    y - FreqSeries object of the coarser series
    index1, index2 - indices of the lowest and highest frequency bins of x \
    that overlap with y (0<=index1<=len(x);1<=index2<=len(x)+1)
    frac1,frac2 - fractional contributions from the above frequency bins
    """
    Nx = len(x.data)
    flowx = x.flow
    deltafx = x.deltaf

    if ( (deltafx <= 0) | (deltafy <= 0) | (Ny <= 0) ):
        raise ValueError, 'bad input argument'
    if ( deltafy < deltafx ):
        raise ValueError, 'deltaf coarse-grain < deltaf fine-grain'
    if ( (flowy - 0.5*deltafy) < (flowx - 0.5*deltafx) ):
        raise ValueError, 'desired coarse-grained start frequency is too low'
    fhighx = flowx + (Nx-1)*deltafx
    fhighy = flowy + (Ny-1)*deltafy
    if ( (fhighy + 0.5*deltafy) > (fhighx + 0.5*deltafx) ):
        raise ValueError, 'desired coarse-rained stop frequency is too high'
    
    i = arange(Ny)
    jlow = int(
        1 + floor((flowy + (i-0.5)*deltafy - flowx - 0.5*deltafx)/deltafx))
    jhigh = int(
        1 + floor((flowy + (i+0.5)*deltafy - flowx - 0.5*deltafx)/deltafx))
    index1 = jlow[0]
    index2 = jhigh[-1]
    fraclow = (flowx + (jlow+0.5)*deltafx - flowy - (i-0.5)*deltafy)/deltafx
    frachigh = (flowy + (i+0.5)*deltafy - flowx - (jhigh-0.5)*deltafx)/deltafx
    frac1 = fraclow[0]
    frac2 = frachigh[-1]

    jtemp = jlow + 1

    midsum = sumTerms(x.data, jtemp, jhigh)
    ya = (deltafx/deltafy)*(x.data[jlow[:-1]]*fraclow[:-1] + 
                            x.data[jhigh[:-1]]*frachigh[:-1] + 
                            midsum[:-1])
    
    if (jhigh[-1] > Nx-1):
        yb = (deltafx/deltafy)*(x.data[jlow[-1]]*fraclow[-1] + midsum[-1])
    else:
        yb = (deltafx/deltafy)*(x.data[jlow[-1]]*fraclow[-1] + 
                                x.data[jhigh[-1]]*frachigh[-1] + 
                                midsum[-1])

    y = FreqSeries()
    y.data = array( list(ya) + [yb] )
    y.flow = flowy
    y.deltaf = deltafy
    return y,index1,index2,frac1,frac2


def windowAndFFT(x,window,fftlength):
    """
    At the moment, this is not equipped to deal with, i think, the heterodyne case properly, because haven't dealt with the issue of NaN yet.  But it should work for the other case.  
    INTPUT:
    x - an instance of TimeSeries
    window - an numpy array of same length as x.data
    fftlength - length of fft
    OUTPUT:
    y - an instance of FreqSeries.  Really, it is, as the name of the function suggests, the result of windowing the time-series x, and Fourier-transforming it.  
    """
    N = len(x.data)
    deltat = x.deltat

    if ( len(x.data) != len(window) ):
        raise Exception, 'size mismatch, length of data ~= that of window'
    if (fftlength < N):
        raise ValueError, 'fftlength < length of data'

    z = zeros(fftlength-N)
    xbar = array( list(x.data*window) + list(z) )
    y_temp = fft(xbar)

    y = FreqSeries()

    if hasattr(x,'fbase'):
        offset = arange(y_temp) - len(y_temp)
        offset = fftshift(offset)
        y_temp = fftshift(y_temp)
        if hasattr(x,'phase') & (x.phase!=0):
            y_temp = y_temp*exp(1j*x.phase)
        y.symmetry = 0
    else:
        y_temp = y_temp[0:floor(fftlength/2+1)] 
        y.symmetry = 1

    y_temp = deltat*y_temp
    y.data = y_temp
    y.deltaf = 1./(deltat*fftlength)
    if hasattr(x,'fbase'):
        y.flow = x.fbase + offset[0]*y.deltaf
    else:
        y.flow = 0

    return y

#some functions for putting together the 'dirty map' and 'covariance matrix'.

def CrossSpectraData(A,E):
    """
    A and E are instances of TimeSeries
    """
    datalength = len(A.data)
    window = hanning(datalength)
    fftlength = datalength
    dataduration = datalength*A.deltat
    Abar = windowAndFFT(A,window,fftlength)
    Ebar = windowAndFFT(E,window,fftlength)

    CAE = FreqSeries()
    CAE.flow = Abar.flow
    CAE.deltaf = Abar.deltaf
    CAE.data = conj(Abar.data)*Ebar.data*2/dataduration
    return CAE

def tsXML_to_psdXML(inxmlpath,outxmldir):
    '''
    Not too nice function for opening a TDIsignal xml file, reading the time-series' inside, calculate the PSDs for each, and then saving them in a similar xml file, but with the time-series replaced by PSDs.  Will have to do for now.
    '''
    inxml = lisaxml.readXML(inxmlpath)
    source = inxml.SourceData[0]
    lisa = inxml.LISAData
    tdiobs = inxml.TDIData[0]
    inxml.close()

    hr = 60.0**2
    patches = int(tdiobs.TimeSeries.Duration/hr)
    specAA = synthlisa.spect(tdiobs.A,tdiobs.TimeSeries.Cadence,patches)
    specEE = synthlisa.spect(tdiobs.E,tdiobs.TimeSeries.Cadence,patches)
    specTT = synthlisa.spect(tdiobs.T,tdiobs.TimeSeries.Cadence,patches)
    f = specAA[:,0]
    pAA = specAA[:,1]
    pEE = specEE[:,1]
    pTT = specTT[:,1]

    pobs = lisaxml.Observable('f,pAA,pEE,pTT',DataType = 'FractionalFrequency')
    pobs.TimeSeries = lisaxml.TimeSeries([f,pAA,pEE,pTT],'f,pAA,pEE,pTT')
    pobs.TimeSeries.Cadence = f[1] - f[0]
    pobs.TimeSeries.Cadence_Unit = 'Hertz'
    pobs.TimeSeries.TimeOffset = 0
    pobs.TimeSeries.TimeOffset_Unit = 'Hertz'

    inxmlname = os.path.basename(inxmlpath)
    outxmlpath = outxmldir + re.sub('\.xml$','',inxmlname) + '-psd.xml'
    outxml = lisaxml.lisaXML(outxmlpath,author='J.Yu')
    outxml.TDIData(pobs)
    outxml.LISAData(lisa)
    outxml.SourceData(source)
    outxml.close()
    return


