#from synthlisa import *
#import matplotlib
#matplotlib.use('Agg')
#from pylab import *
import synthlisa
from spharm import Spharmt, getspecindx
import numpy as np
import scipy as sp
import scipy.linalg as sp_linalg


def get_tdi_lnklist( tdi_type , tdi_generation , conv_ref ) :
    """
    INPUT:
    tdi_type --- 'Michelson' , 'Sagnac' , 'Beacon'
    tdi_generation --- 'G0' , 'G1' , 'Gm' , 'G2'
    conv_ref --- conventional reference spacecraft: '1', '2' or '3'
    OUTPUT:
    tdi_lnklist --- tuple containing info about links in the tdi observable
    """
    if tdi_type not in [ 'Michelson' , 'Sagnac' ] :
        raise InputError , 'TDI type not known!'
    if tdi_generation not in [ 'G0' , 'G1' , 'Gm' , 'G2' ] :
        raise InputError , 'TDI generation not known!'
    if conv_ref not in [ '1' , '2' , '3' ] :
        raise InputError , "Reference spacecraft must be either '1', '2', or '3'!"
    
    if tdi_type == 'Michelson' :
        if tdi_generation == 'G0' :
            if conv_ref == '1' :
                tdi_lnklist = ({'lnk':'123','sgn':-1,'tdi':[(-1,-2)]},
                               {'lnk':'3-21','sgn':-1,'tdi':[]},
                               {'lnk':'231','sgn':+1 ,'tdi':[]} , 
                               {'lnk':'1-32','sgn':+1,'tdi':[(-1,3)]}) 
            elif conv_ref == '2' :
                tdi_lnklist = ({'lnk':'231','sgn':-1,'tdi':[(-1,-3)]},
                               {'lnk':'1-32','sgn':-1,'tdi':[]},
                               {'lnk':'312','sgn':+1,'tdi':[]}, 
                               {'lnk':'2-13','sgn':+1,'tdi':[(-1,1)]})
            elif conv_ref == '3' :
                tdi_lnklist = ({'lnk':'312','sgn':-1,'tdi':[(-1,-1)]},
                               {'lnk':'2-13','sgn':-1,'tdi':[]},
                               {'lnk':'123','sgn':+1,'tdi':[]}, 
                               {'lnk':'3-21','sgn':+1,'tdi':[(-1,2)]})
        elif tdi_generation == 'Gm' :
            if conv_ref == '1' :
                tdi_lnklist = ({'lnk':'1-32','sgn':1,'tdi':[(-1,-2),(-1,2),(-1,3)]},
                               {'lnk':'231','sgn':1,'tdi':[(-1,-2),(-1,2)]},
                               {'lnk':'123','sgn':1,'tdi':[(-1,-2)]},
                               {'lnk':'3-21','sgn':1,'tdi':[]},
                               {'lnk':'231','sgn':-1,'tdi':[]}, 
                               {'lnk':'1-32','sgn':-1,'tdi':[(-1,3)]},
                               {'lnk':'3-21','sgn':-1,'tdi':[(-1,3),(-1,-3)]},
                               {'lnk':'123','sgn':-1,'tdi':[(-1,3),(-1,-3),(-1,-2)]}) 
            elif conv_ref == '2' :
                tdi_lnklist = ({'lnk':'2-13','sgn':1,'tdi':[(-1,-3),(-1,3),(-1,1)]},
                               {'lnk':'312','sgn':1,'tdi':[(-1,-3),(-1,3)]} ,
                               {'lnk':'231','sgn':1,'tdi':[(-1,-3)]} ,
                               {'lnk':'1-32','sgn':1,'tdi':[]},
                               {'lnk':'312','sgn':-1,'tdi':[]}, 
                               {'lnk':'2-13','sgn':-1,'tdi':[(-1,1)]} ,
                               {'lnk':'1-32','sgn':-1,'tdi':[(-1,1),(-1,-1)]} ,
                               {'lnk':'231','sgn':-1,'tdi':[(-1,1),(-1,-1),(-1,-3)]}) 
            elif conv_ref == '3' :
                tdi_lnklist = ({'lnk':'3-21','sgn':1,'tdi':[(-1,-1),(-1,1),(-1,2)]},
                               {'lnk':'123','sgn':1,'tdi':[(-1,-1),(-1,1)]},
                               {'lnk':'312','sgn':1,'tdi':[(-1,-1)]},
                               {'lnk':'2-13','sgn':1,'tdi':[]},
                               {'lnk':'123','sgn':-1,'tdi':[]}, 
                               {'lnk':'3-21','sgn':-1,'tdi':[(-1,2)]},
                               {'lnk':'2-13','sgn':-1,'tdi':[(-1,2),(-1,-2)]},
                               {'lnk':'312','sgn':-1,'tdi':[(-1,2),(-1,-2),(-1,-1)]}) 
        elif tdi_generation == 'G2' :
            if conv_ref == '1' :
                tdi_lnklist = ({'lnk':'1-32','sgn':1,'tdi':[(-1,3),(-1,-3),(-1,-2),(-1,2),(-1,-2),(-1,2),(-1,3)]},
                               {'lnk':'231','sgn':1,'tdi':[(-1,3),(-1,-3),(-1,-2),(-1,2),(-1,-2),(-1,2)]},
                               {'lnk':'123','sgn':1,'tdi':[(-1,3),(-1,-3),(-1,-2),(-1,2),(-1,-2)]},
                               {'lnk':'3-21','sgn':1,'tdi':[(-1,3),(-1,-3),(-1,-2),(-1,2)]},
                               {'lnk':'123','sgn':1,'tdi':[(-1,3),(-1,-3),(-1,-2)]},
                               {'lnk':'3-21','sgn':1,'tdi':[(-1,3),(-1,-3)]},
                               {'lnk':'1-32','sgn':1,'tdi':[(-1,3)]},
                               {'lnk':'231','sgn':1,'tdi':[]},
                               {'lnk':'3-21','sgn':-1,'tdi':[]},
                               {'lnk':'123','sgn':-1,'tdi':[(-1,-2)]},
                               {'lnk':'231','sgn':-1,'tdi':[(-1,-2),(-1,2)]},
                               {'lnk':'1-32','sgn':-1,'tdi':[(-1,-2),(-1,2),(-1,3)]},
                               {'lnk':'231','sgn':-1,'tdi':[(-1,-2),(-1,2),(-1,3),(-1,-3)]},
                               {'lnk':'1-32','sgn':-1,'tdi':[(-1,-2),(-1,2),(-1,3),(-1,-3),(-1,3)]},
                               {'lnk':'3-21','sgn':-1,'tdi':[(-1,-2),(-1,2),(-1,3),(-1,-3),(-1,3),(-1,-3)]},
                               {'lnk':'123','sgn':-1,'tdi':[(-1,-2),(-1,2),(-1,3),(-1,-3),(-1,3),(-1,-3),(-1,-2)]})
            elif conv_ref == '2' :
                tdi_lnklist = ({'lnk':'2-13','sgn':1,'tdi':[(-1,1),(-1,-1),(-1,-3),(-1,3),(-1,-3),(-1,3),(-1,1)]},
                               {'lnk':'312','sgn':1,'tdi':[(-1,1),(-1,-1),(-1,-3),(-1,3),(-1,-3),(-1,3)]},
                               {'lnk':'231','sgn':1,'tdi':[(-1,1),(-1,-1),(-1,-3),(-1,3),(-1,-3)]},
                               {'lnk':'1-32','sgn':1,'tdi':[(-1,1),(-1,-1),(-1,-3),(-1,3)]},
                               {'lnk':'231','sgn':1,'tdi':[(-1,1),(-1,-1),(-1,-3)]},
                               {'lnk':'1-32','sgn':1,'tdi':[(-1,1),(-1,-1)]},
                               {'lnk':'2-13','sgn':1,'tdi':[(-1,1)]},
                               {'lnk':'312','sgn':1,'tdi':[]} ,
                               {'lnk':'1-32','sgn':-1,'tdi':[]} ,
                               {'lnk':'231','sgn':-1,'tdi':[(-1,-3)]} ,
                               {'lnk':'312','sgn':-1,'tdi':[(-1,-3),(-1,3)]},
                               {'lnk':'2-13','sgn':-1,'tdi':[(-1,-3),(-1,3),(-1,1)]},
                               {'lnk':'312','sgn':-1,'tdi':[(-1,-3),(-1,3),(-1,1),(-1,-1)]},
                               {'lnk':'2-13','sgn':-1,'tdi':[(-1,-3),(-1,3),(-1,1),(-1,-1),(-1,1)]},
                               {'lnk':'1-32','sgn':-1,'tdi':[(-1,-3),(-1,3),(-1,1),(-1,-1),(-1,1),(-1,-1)]},
                               {'lnk':'231','sgn':-1,'tdi':[(-1,-3),(-1,3),(-1,1),(-1,-1),(-1,1),(-1,-1),(-1,-3)]})
            elif conv_ref == '3' :
                tdi_lnklist = ({'lnk':'3-21','sgn':1,'tdi':[(-1,2),(-1,-2),(-1,-1),(-1,1),(-1,-1),(-1,1),(-1,2)]},
                               {'lnk':'123','sgn':1,'tdi':[(-1,2),(-1,-2),(-1,-1),(-1,1),(-1,-1),(-1,1)]},
                               {'lnk':'312','sgn':1,'tdi':[(-1,2),(-1,-2),(-1,-1),(-1,1),(-1,-1)]},
                               {'lnk':'2-13','sgn':1,'tdi':[(-1,2),(-1,-2),(-1,-1),(-1,1)]},
                               {'lnk':'312','sgn':1,'tdi':[(-1,2),(-1,-2),(-1,-1)]},
                               {'lnk':'2-13','sgn':1,'tdi':[(-1,2),(-1,-2)]},
                               {'lnk':'3-21','sgn':1,'tdi':[(-1,2)]},
                               {'lnk':'123','sgn':1,'tdi':[]},
                               {'lnk':'2-13','sgn':-1,'tdi':[]},
                               {'lnk':'312','sgn':-1,'tdi':[(-1,-1)]},
                               {'lnk':'123','sgn':-1,'tdi':[(-1,-1),(-1,1)]},
                               {'lnk':'3-21','sgn':-1,'tdi':[(-1,-1),(-1,1),(-1,2)]},
                               {'lnk':'123','sgn':-1,'tdi':[(-1,-1),(-1,1),(-1,2),(-1,-2)]},
                               {'lnk':'3-21','sgn':-1,'tdi':[(-1,-1),(-1,1),(-1,2),(-1,-2),(-1,2)]},
                               {'lnk':'2-13','sgn':-1,'tdi':[(-1,-1),(-1,1),(-1,2),(-1,-2),(-1,2),(-1,-2)]},
                               {'lnk':'312','sgn':-1,'tdi':[(-1,-1),(-1,1),(-1,2),(-1,-2),(-1,2),(-1,-2),(-1,-1) ] } )
    elif tdi_type == 'Sagnac' :
        if tdi_generation == 'G1' :
            if conv_ref == '1' :
                tdi_lnklist = ( {'lnk':'123','sgn':-1,'tdi':[(-1,3),(-1,1)]},
                                {'lnk':'312','sgn':-1,'tdi':[(-1,3)]},
                                {'lnk':'231','sgn':-1,'tdi':[]},
                                {'lnk':'3-21','sgn':+1,'tdi':[]},
                                {'lnk':'2-13','sgn':+1,'tdi':[(-1,-2)]},
                                {'lnk':'1-32','sgn':+1,'tdi':[(-1,-2),(-1,-1)]} )
            elif conv_ref == '2' :
                tdi_lnklist = ( {'lnk':'231','sgn':-1,'tdi':[(-1,1),(-1,2)]},
                                {'lnk':'123','sgn':-1,'tdi':[(-1,1)]},
                                {'lnk':'312','sgn':-1,'tdi':[]},
                                {'lnk':'1-32','sgn':+1,'tdi':[]},
                                {'lnk':'3-21','sgn':+1,'tdi':[(-1,-3)]},
                                {'lnk':'2-13','sgn':+1,'tdi':[(-1,-3),(-1,-2)]} )
            elif conv_ref == '3' :
                tdi_lnklist = ( {'lnk':'312','sgn':-1,'tdi':[(-1,2),(-1,3)]},
                                {'lnk':'231','sgn':-1,'tdi':[(-1,2)]},
                                {'lnk':'123','sgn':-1,'tdi':[]},
                                {'lnk':'2-13','sgn':+1,'tdi':[]},
                                {'lnk':'1-32','sgn':+1,'tdi':[(-1,-1)]},
                                {'lnk':'3-21','sgn':+1,'tdi':[(-1,-1),(-1,-3)]} )
        elif tdi_generation == 'Gm' :
            if conv_ref == '1' :
                tdi_lnklist = ({'lnk':'123','sgn':-1,'tdi':[(-1,-2),(-1,-1),(-1,-3),(-1,3),(-1,1)]},
                               {'lnk':'312','sgn':-1,'tdi':[(-1,-2),(-1,-1),(-1,-3),(-1,3)]},
                               {'lnk':'231','sgn':-1,'tdi':[(-1,-2),(-1,-1),(-1,-3)]},
                               {'lnk':'1-32','sgn':-1,'tdi':[(-1,-2),(-1,-1)]},
                               {'lnk':'2-13','sgn':-1,'tdi':[(-1,-2)]},
                               {'lnk':'3-21','sgn':-1,'tdi':[]},
                               {'lnk':'231','sgn':+1,'tdi':[]},
                               {'lnk':'312','sgn':+1,'tdi':[(-1,3)]},
                               {'lnk':'123','sgn':+1,'tdi':[(-1,3),(-1,1)]},
                               {'lnk':'3-21','sgn':+1,'tdi':[(-1,3),(-1,1),(-1,2)]},
                               {'lnk':'2-13','sgn':+1,'tdi':[(-1,3),(-1,1),(-1,2),(-1,-2)]},
                               {'lnk':'1-32','sgn':+1,'tdi':[(-1,3),(-1,1),(-1,2),(-1,-2),(-1,-1)]})
            elif conv_ref == '2' :
                tdi_lnklist = ({'lnk':'231','sgn':-1,'tdi':[(-1,-3),(-1,-2),(-1,-1),(-1,1),(-1,2)]},
                               {'lnk':'123','sgn':-1,'tdi':[(-1,-3),(-1,-2),(-1,-1),(-1,1)]},
                               {'lnk':'312','sgn':-1,'tdi':[(-1,-3),(-1,-2),(-1,-1)]},
                               {'lnk':'2-13','sgn':-1,'tdi':[(-1,-3),(-1,-2)]},
                               {'lnk':'3-21','sgn':-1,'tdi':[(-1,-3)]},
                               {'lnk':'1-32','sgn':-1,'tdi':[]},
                               {'lnk':'312','sgn':+1,'tdi':[]},
                               {'lnk':'123','sgn':+1,'tdi':[(-1,1)]},
                               {'lnk':'231','sgn':+1,'tdi':[(-1,1),(-1,2)]},
                               {'lnk':'1-32','sgn':+1,'tdi':[(-1,1),(-1,2),(-1,3)]},
                               {'lnk':'3-21','sgn':+1,'tdi':[(-1,1),(-1,2),(-1,3),(-1,-3)]},
                               {'lnk':'2-13','sgn':+1,'tdi':[(-1,1),(-1,2),(-1,3),(-1,-3),(-1,-2)]})
            elif conv_ref == '3' :
                tdi_lnklist = ({'lnk':'312','sgn':-1,'tdi':[(-1,-1),(-1,-3),(-1,-2),(-1,2),(-1,3)]},
                               {'lnk':'231','sgn':-1,'tdi':[(-1,-1),(-1,-3),(-1,-2),(-1,2)]},
                               {'lnk':'123','sgn':-1,'tdi':[(-1,-1),(-1,-3),(-1,-2)]},
                               {'lnk':'3-21','sgn':-1,'tdi':[(-1,-1),(-1,-3)]},
                               {'lnk':'1-32','sgn':-1,'tdi':[(-1,-1)]},
                               {'lnk':'2-13','sgn':-1,'tdi':[]},
                               {'lnk':'123','sgn':+1,'tdi':[]},
                               {'lnk':'231','sgn':+1,'tdi':[(-1,2)]},
                               {'lnk':'312','sgn':+1,'tdi':[(-1,2),(-1,3)]},
                               {'lnk':'2-13','sgn':+1,'tdi':[(-1,2),(-1,3),(-1,1)]},
                               {'lnk':'1-32','sgn':+1,'tdi':[(-1,2),(-1,3),(-1,1),(-1,-1)]},
                               {'lnk':'3-21','sgn':+1,'tdi':[(-1,2),(-1,3),(-1,1),(-1,-1),(-1,-3)]})
    return tdi_lnklist




class planeGWresponse(object):
    def __init__(self,LISA,elat,elon,pol):
        """
        INPUT:
        LISA --- synthetic LISA object
        elat --- ecliptic latitude [radians]
        elon --- ecliptic longitude [radians]
        pol --- polarisation angle [radians]
        """
        self.Lisa = LISA
        self.wave = synthlisa.SimpleMonochromatic(1e-3,0,1,2,elat,elon,pol)
        self.elat = elat
        self.elon = elon
        self.pol = pol
        return

    def getk(self):
        """
        returns components of the unit vector in the direction of propagation
        OUPUT:
        numpy array containing the above
        """
        return np.array(self.wave.putk())

    def getepol(self,polarisation='plus'):
        """
        INPUT:
        polarisation --- 'plus' or 'cross' GW polarisation
        OUPUT:
        numpy array --- components of GW polarisation basis tensor
        """
        try:
            if polarisation == 'plus':
                return np.array(self.wave.putep(self.elat,self.elon,self.pol))
            elif polarisation == 'cross':
                return np.array(self.wave.putec(self.elat,self.elon,self.pol))
            else:
                raise ValueError
        except ValueError:
            print "polarisation must be either 'plus' or 'cross'"
        return

    def getp(self,n,t=0):
        """
        INPUT:
        n --- spacecraft number 1 , 2 or 3
        t --- time [s]
        OUTPUT:
        numpy array --- components of the position vector of spacecraft n at time t
        """
        return np.array(self.Lisa.putp(n,t)) 
    def getn(self,n,t=0):
        """
        INPUT:
        n --- link number 1, 2, 3, -1, -2, or -3
        t --- time [s], reception time
        OUTPUT:
        numpy array --- components of the unit vector along link n at time t
        """
        return np.array(self.Lisa.putn(n,t)) 

    def armlength(self,n,t=0):
        """
        INPUT:
        n --- link number 1, 2, 3, -1, -2, or -3
        t --- time [s] , reception time
        OUTPUT:
        optical path length along link n at time t. [s]
        """
        return self.Lisa.armlength(n,t)

    def getfn(self,n,t=0):
        """
        INPUT:
        n --- link number 1, 2, 3, -1, -2 or -3
        t -- time [s] , reception time
        OUTPUT:
        characteristic transfer frequency of link n at time t [Hz]
        """
        return 1/(2*np.pi*self.Lisa.armlength(n,t))

    def getlink_arm(self,link):
        """
        INPUT:
        link --- string indicating one of the following links:
        '123','312','231','1-32','2-13','3-21'
        OUTPUT:
        link number, i.e the number in the middle
        """
        link_dictionary = {'123':2,'312':1,'231':3,
                           '1-32':-3,'2-13':-1,'3-21':-2}
        return link_dictionary[link]

    def getlink_receiver(self,link):
        """
        INPUT:
        link --- string indicating one of:  '123','312','231','1-32','2-13','3-21'
        OUTPUT:
        receiving spacecraft number, i.e the last number
        """
        return int(link[-1])
        
    def linkF(self,link,f,t=0,polarisation='plus') :
        """
        The response function of a one-way Doppler measurement.
        The event of reception is taken to be the reference.
        INPUT:
        f --- frequency/frequencies
        link --- '123','312','231','1-32','2-13','3-21'
        t --- time [s] 
        polarisation --- 'plus' or 'cross' polarisation
        OUTPUT:
        response function of the link
        """

        l = self.getlink_arm(link) 
        n = self.getn(l,t) 
        L = self.armlength(l,t)
        fl = self.getfn(l,t)
        k = self.getk()
        E = self.getepol(polarisation)
        kn = np.dot(k,n)
        nEn = np.dot(n,np.dot(E,n))
        argu = (1-kn)*f/fl/2
        Famp = -1j*(f/fl/2)*nEn*np.sinc(argu/np.pi) #fractional frequency
#        Famp = -1j*(f/fl/2)*nEn*sinc(argu/pi) / (1j*2*pi*f) #phase difference       
        Fphse = -argu
        F = Famp*np.exp(1j*Fphse)
        return F


    def timedelay( self , tdi , t=0 ) :
        """
        INPUT:
        tdi --- list containing a sequence of time-delay/advance tuples [ ( +/-1 , Ll ) , ... ]
                the delays and advances are counted in the same order as the list tdi
        t --- time from which to start counting the sequence of delays and advances
        OUTPUT:
        net elapse of time from t
        """
        netelapse = 0
        for delay in tdi:
            elapse = delay[0] * self.armlength( delay[1] , t + netelapse )
            netelapse += elapse
        return netelapse

    def tdiF( self , tdi_type , tdi_generation , conv_ref ,
              f , t=0 , polarisation='plus' , ref='1' ) :
        """
        INPUT:
        tdi_type --- type of TDI observable: 'Michelson' , 'Sagnac' etc.
        tdi_generation --- TDI generation: 'G0','G1','Gm','G2'
        conv_ref --- conventional reference spacecraft number: X <-> '1' , Y <-> '2' etc.
        f --- frequency/frequencies [Hz]
        t --- time [s]
              This is the time of reception of the links with ZERO time-delay in the TDI observable as defined in tdi_lnklist()
        polarisation --- GW polarisation ( default is + )
        ref --- user-chosen reference spacecraft number (this does not define the TDI).
                The position of this spacecraft, at time t, is used as the reference position.
        OUTPUT:
        tdiF --- GW response function of TDI observable
        """
        tdi_lnklist = get_tdi_lnklist( tdi_type , tdi_generation , conv_ref )
        refp = self.getp( int( ref ) , t )
        F = np.zeros( len(f) , dtype = complex ) ###NOTE complex64 straight works here because zeros and not pylab.zeros 
        for lnk in tdi_lnklist:
            td = self.timedelay( lnk['tdi'] , t )
            tlnk = t + td
            lnkp = self.getp( int(lnk['lnk'][-1]) , tlnk )
            kp = np.dot( self.getk() , lnkp-refp )
            lnkF = self.linkF( lnk['lnk'] , f , tlnk , polarisation )
            F += lnk['sgn'] * lnkF * np.exp( 1j*2*np.pi*f*(td-kp) )
        tdiF = F
        return tdiF

    def optimal_tdiF( self , tdi_type , tdi_generation , whichtdi ,
                      f , t=0 , polarisation='plus' , ref='1' ) :
        """
        INPUT:
        tdi_type --- type of TDI used as basis for optimal TDI observables
        tdi_generation -- TDI generation: 'G0','G1','Gm', or 'G2'
        whichtdi --- name of optimal TDI observable: 'A', 'E' or 'T'
        f --- frequency/frequencies
        t --- time [s]
              This is the time of reception of the links with ZERO time-delay in the TDI observable as defined in tdi_lnklist()
        polarisation --- GW polarisation ( default is + )
        ref --- user-chosen reference spacecraft number
                The position of this spacecraft, at time t, is used as the reference position.
        OUTPUT:
        tdiF --- GW response function of optimal TDI observable
        """
        if tdi_type not in [ 'Michelson' , 'Sagnac' ] :
            raise InputError , 'Unrecognised TDI type as the basis for the optimal observables!'
        if whichtdi not in [ 'A' , 'E' , 'T' ] :
            raise InputError , 'Unrecognised optimal TDI name!'
        tdiFs = [ self.tdiF( tdi_type , tdi_generation , i , f , t , polarisation , ref )
                  for i in ['1','2','3'] ]        
        if tdi_type == 'Sagnac' :
            alphaF , betaF , gammaF = tdiFs
            if whichtdi == 'A' :
                tdiF = ( gammaF - alphaF ) / np.sqrt(2)
            elif whichtdi == 'E' :
                tdiF = ( alphaF - 2*betaF + gammaF ) / np.sqrt(6)
            elif whichtdi == 'T' :
                tdiF = ( alphaF + betaF + gammaF ) / np.sqrt(3)
        elif tdi_type == 'Michelson' :
            XF , YF , ZF = tdiFs
            if whichtdi == 'A' :
                tdiF = ( 2*XF - YF - ZF ) / 3
            elif whichtdi == 'E' :
                tdiF = ( ZF - YF )/ np.sqrt(3)
            elif whichtdi == 'T' :
                tdiF = ( XF + YF + ZF  ) / 3
        return tdiF

    def tdiORF( self , tdiI , tdiJ , f , t=0 ) :
        """
        INPUT:
        tdiI --- tuple (tdi_type,tdi_generation,con_ref/whichtdi,ref) for TDI observ. I
        tdiJ --- tuple (tdi_type,tdi_generation,con_ref/whichtdi,ref) for TDI observ. J
        f --- frequency/frequencies [Hz]
        t --- time [s]
        OUTPUT:
        orf --- overlap reduction function between tdiI and tdiJ
        """
        ps , tdiFs = [] , []
        for k , tdi in enumerate( [ tdiI , tdiJ ] ) :
            tdi_type , tdi_generation , whichtdi , ref = tdi
            ps  += [ self.getp( int( ref ) , t ) ]
            if whichtdi in [ '1' , '2' , '3' ] :
                tdiFs += [ [ self.tdiF( tdi_type , tdi_generation , whichtdi ,
                                        f , t , pol , ref )
                             for pol in ['plus','cross'] ] ]
            elif whichtdi in [ 'A' , 'E' , 'T' ] :
                tdiFs += [ [ self.optimal_tdiF( tdi_type , tdi_generation , whichtdi ,
                                                f , t , pol , ref )
                           for pol in ['plus','cross'] ] ]
            else:
                raise InputError , '%d th whichtdi is not recognised!' % (k+1)
        kp = np.dot( self.getk() , ps[0] - ps[1] )
        orf = ( np.conj( tdiFs[0][0] )*tdiFs[1][0] +
                np.conj( tdiFs[0][1] )*tdiFs[1][1] ) * np.exp( 1j * 2*np.pi * f * kp )
        return orf 



class mySpharmt(Spharmt):
    def __init__(self,nlon,nlat,rsphere=1e2,gridtype='regular',legfunc='stored'):
        Spharmt.__init__(self,nlon,nlat,rsphere,gridtype,legfunc)
        if self.gridtype=='gaussian':
            raise Exception, "Can't deal with Gaussian grid at the moment."
        else:
            deltalat = 180./(nlat-1)
            self.lats = 90 - deltalat*np.arange(nlat)
        deltalon = 360./nlon
        self.lons = deltalon*np.arange(nlon)
        return



class LISA_in_the_Sky(object):
    def __init__(self,LISA,mySpharmt):
        self.Lisa = LISA
        self.sky = mySpharmt
        return

    def get_pixels( self , whichfunc , *args , **kargs ) :
        """
        INPUT:
        whichfunc --- which function to evaluate over the sky: 'linkF' , 'tdiF' , 'optimal_tdiF' or 'tdiORF'
        args/kargs --- these are the same as the corresponding methods in planeGWresponse()
        """
        nlat = self.sky.nlat ; lats = self.sky.lats
        nlon = self.sky.nlon ; lons = self.sky.lons
        firstpixel = True
        for x , elatd in enumerate(lats):
            for y , elond in enumerate(lons):
                elat , elon = elatd * np.pi / 180 , elond * np.pi / 180
#                print 'myLISAmodule.py: LISA_in_the_Sky.get_pixels() working on ( elat , elon ) = ( %f , %f )' % ( np.degrees( elat ) , np.degrees( elon ) )
                response = planeGWresponse( self.Lisa , elat , elon , 0 )
                if whichfunc == 'linkF' :
                    if len( args ) < 2 :
                        raise InputError, 'You must at least specify the slr link and the frequencies!'
                    if firstpixel :
                        func = np.zeros( ( nlat , nlon , len( args[1] ) ) , dtype = complex ) ; firstpixel = False
                    if len( args ) == 4 :
                        func[x,y,:] = response.linkF( *args )
                    elif len( args ) == 3 :
                        func[x,y,:] = response.linkF( polarisation=kargs['polarisation'] , *args )
                    elif len( args ) == 2 :
                        func[x,y,:] = response.linkF( polarisation=kargs['polarisation'] , t=kargs['t'] , *args )
                elif whichfunc == 'tdiF' :
                    if len( args ) < 4 :
                        raise InputError, 'You must at least specify tditype, tdigen, whichtdi and f!'
                    if firstpixel :
                        func = np.zeros( ( nlat , nlon , len( args[3] ) ) , dtype = complex ) ; firstpixel = False
                    if len( args ) == 7 :
                        func[ x , y , : ] = response.tdiF( *args )
                    elif len( args ) == 6 :
                        func[ x , y , : ] = response.tdiF( ref = kargs['ref'] , *args )
                    elif len( args ) == 5 :
                        func[ x , y , : ] = response.tdiF( polaristaion = kargs['polarisation'] , ref = kargs['ref'] , *args )
                    elif len( args ) == 4 :
                        func[ x , y , : ] = response.tdiF( t = kargs['t'] , polarisation = kargs['polarisation'] , ref = kargs['ref'] , *args )
                elif whichfunc == 'optimal_tdiF' :
                    if len( args ) < 4 :
                        raise InputError, 'You must at least specify tditype, tdigen, whichtdi and f!'
                    if firstpixel :
                        func = np.zeros( ( nlat , nlon , len( args[3] ) ) , dtype = complex ) ; firstpixel = False
                    if len( args ) == 7 :
                        func[ x , y , : ] = response.optimal_tdiF( *args )
                    elif len( args ) == 6 :
                        func[ x , y , : ] = response.optimal_tdiF( ref = kargs['ref'] , *args )
                    elif len( args ) == 5 :
                        func[ x , y , : ] = response.optimal_tdiF( polaristaion = kargs['polarisation'] , ref = kargs['ref'] , *args )
                    elif len( args ) == 4 :
                        func[ x , y , : ] = response.optimal_tdiF( t = kargs['t'] , polarisation = kargs['polarisation'] , ref = kargs['ref'] , *args)
                elif whichfunc == 'tdiORF' :
                    if len( args ) < 3 :
                        raise InputError , 'You must at least specify tdiI, tdiJ and f!'
                    if firstpixel :
                        func = np.zeros( ( nlat , nlon , len( args[2] ) ) , dtype = complex ) ; firstpixel = False
                    if len( args ) == 4 :
                        func[ x , y , : ] = response.tdiORF( *args )
                    elif len( args ) == 3 :
                        func[ x , y , : ] = response.tdiORF( t = kargs['t'] , *args )
                else:
                    raise InputError , 'Unrecognised function %s' % func
        return func 

    def get_SpHs( self , ntrunc , whichfunc , *args , **kargs ) :
        """
        INPUT:
        whichfunc --- which function: 'linkF' , 'tdiF' , 'optimal_tdiF' , or 'tdiORF'
        ntrunc --- maximum l up to which to expand the angular function in terms of SpH
        OUTPUT:
        specr --- multipole moments of the real part of orf for positive m
        speci --- multipole moments of the imag part of orf for positive m        
        """
        sky = self.sky
        func = self.get_pixels( whichfunc , *args , **kargs )
        nlm = ( ntrunc + 1 )*( ntrunc + 2 ) / 2
        numf = np.size( func , 2 )
        specr = np.zeros( ( nlm , numf ) , dtype = complex )
        speci = np.zeros( ( nlm , numf ) , dtype = complex )
        for i in range( numf ) :
            print "myLISAmodule.py: lisa_in_the_Sky.getSpHs() working on %d th frequency of %d" % ( i+1 , numf )
            specr[ : , i ] = sky.grdtospec( np.real( func[:,:,i]) , ntrunc )
            speci[ : , i ] = sky.grdtospec( np.imag( func[:,:,i]) , ntrunc )
        self.ntrunc = ntrunc
        return specr , speci        




def getMLvec(ntrunc,m='pn'):
    if ntrunc==0:
        M = [0]
        L = [0]
    elif m=='pn':    
        m,l = getspecindx(ntrunc)
        M = list(m) + list(-m[ntrunc+1:])
        L = list(l) + list(l[ntrunc+1:])
    elif m=='p':
        M,L = getspecindx(ntrunc)
    return zip(M,L)


def ORFfromFiles(realrpath, realipath, imagrpath, imagipath, sign='pn'):
    if not (sign=='p' or sign=='pn'):
        raise ValueError,"Keyword parameter sign must either be 'p' or 'pn'."
    realr = load(realrpath)
    reali = load(realipath)
    imagr = load(imagrpath)
    imagi = load(imagipath)
    try:
        assert shape(realr)==shape(reali)==shape(imagr)==shape(imagi)
    except AssertionError:
        print "The four arrays uploaded need to have the same dimensions."
    ntrunc = int( -1.5 + np.sqrt( 8*shape(realr)[0] + 1 ) / 2  ) 
    numf = shape(realr)[1]
    if (sign == 'p'):
        indxpn = getMLvec(ntrunc,'p')
        p = realr + 1j*reali
        q = imagr + 1j*imagi
        g = np.zeros( ( len(indxpn),numf ) , dtype = complex )
        g = p + 1j*q
    else:
        indxp  = getMLvec(ntrunc,'p')
        indxpn = getMLvec(ntrunc)
        p = realr + 1j*reali
        q = imagr + 1j*imagi
        g = np.zeros( ( len(indxpn),numf ) , dtype = complex )
        for i,ml in enumerate(indxpn):
            m = ml[0]
            l = ml[1]
            if m >= 0:
                k = indxp.index(ml)
                g[i,:] = p[k,:] + 1j*q[k,:]
            else:
                ml = (-m,l)
                k = indxp.index(ml)
                g[i,:] = (-1)**m*( np.conj(p[k,:]) + 1j*np.conj(q[k,:]) )
    return indxpn,g



def get_pmNSD( f ) :
    """
    returns noise spectral density of standard proof-mass noise
    INPUT:
    f --- frequency/frequencies
    OUTPUT:
    psdpm --- noise spectral density [Hz^{-1}]
    """
    psdpm = 2.54e-48 * f**-2
    return psdpm

def get_opNSD( f ) :
    """
    returns noise spectral density of standard optical path noise
    INPUT:
    f --- frequency/frequencies
    OUTPUT:
    psdop --- noise spectral density [Hz^{-1}]
    """
    psdop = 1.76e-37 * f**2
    return psdop

def get_tdiNSD( tditype , tdigen , tdiI , tdiJ , f ) :
    """
    Noise spectral density for TDI observables
    INPUT :
    tdi_type --- type of TDI used as basis for optimal TDI observables
    tdi_generation -- TDI generation: 'G0','G1','Gm', or 'G2'
    tdiI --- which TDI observable: '1' , '2' , '3' , 'A', 'E', or 'T'
    tdiJ --- which TDI observable: '1' , '2' , '3' , 'A', 'E', or 'T'
    f --- frequency/frequencies
    OUTPUT :
    noise spectral density
    """
    if tditype not in [ 'Michelson' , 'Sagnac' ] :
        raise InputError , 'Unrecognised TDI type as the basis for the optimal observables!'
    if tdigen not in [ 'G0' , 'Gm' , 'G1' , 'G2' ] :
        raise InputError , 'Unrecognised TDI generation!'
    if ( tdiI or tdiJ ) not in [ 'A' , 'E' , 'T' , '1' , '2' , '3' ] :
        raise InputError , 'Unrecognised optimal TDI name! '

    L = 16.6782 
    Ppm = get_pmNSD( f )
    Pop = get_opNSD( f )

    if tditype == 'Michelson' :
        if tdigen == 'G0' :
            if ( tdiI , tdiJ ) in [ ('1','1') , ('2','2') , ('3','3') ] :
                nsd = 4 * ( Pop + ( 3 + np.cos( 4*np.pi*L * f ) ) * Ppm )
            elif ( tdiI , tdiJ ) in [ ('1','2') , ('1','3') , ('2','1') , ('2','3') , ('3','1') , ('3','2') ] :
                nsd = - 2 * np.cos( 2*np.pi*L * f ) * ( Pop + 4 * Ppm )
            elif ( tdiI , tdiJ ) in [ ('A','A') , ('E','E') ] :
                nsd = (4./3) * ( ( 2 + np.cos(2*np.pi*L*f) ) * Pop + 2 * ( 3 + 2*np.cos(2*np.pi*L*f) + np.cos(4*np.pi*L*f) ) * Ppm )
            elif ( tdiI , tdiJ ) == ( 'T' , 'T' ) :
                nsd = (8./3) * np.sin( np.pi*L*f )**2 * ( Pop + 2 * ( 1 - np.cos(2*np.pi*L*f) ) * Ppm )
            elif ( tdiI , tdiJ ) in [ ('A','E') , ('A','T') , ('E','A') , ('E','T') , ('T','A') , ('T','E') ] :
                nsd = 0. * f
        if tdigen == 'Gm' :
            nsd = get_tdiNSD( tditype , 'G0' , tdiI , tdiJ , f ) * 4 * np.sin(2*np.pi*L*f)**2 
        if tdigen == 'G2' :
            nsd = get_tdiNSD( tditype , 'G0' , tdiI , tdiJ , f ) * 16 * np.sin(4*np.pi*L*f)**2 * np.sin(2*np.pi*L*f)**2
    return nsd



def get_tdiNoise_freq_domain_1NSD( tditype , tdigen , tdiI , tdiJ ,
                                   stime , duration , inittime , seed ) :
    """
    returns the noise time-series given its power spectral density
    INPUT:
    tditype --- which type of TDI observables
    tdigen --- which generation of TDI
    tdiI --- which TDI observable. 'A' , 'E' , 'T' , '1' , '2' , '3'
    tdiJ --- which TDI observable. 'A' , 'E' , 'T' , '1' , '2' , '3'
    stime --- sample time of the noise time-series
    duration --- duration of the noise time-series
    inittime --- initial time of the noise time-series
    seed --- seed for the noise
    OUTPUT:
    t --- time [s]
    n --- noise time-series
    """
    np.random.seed( seed )
    if tdiI != tdiJ :
        raise InputError , 'It only makes sense to simulate noise from its PSD, so please make sure tdiI and tdiJ are the same!'
    N = np.round( duration / stime )
    t = inittime + stime * np.arange( N )
    fs = 1. / stime ; df = 1 / ( N*stime )
    if N % 2 == 0 :
        Nf = N/2 - 1
    else :
        Nf = ( N - 1 ) / 2
    f = df * np.arange( 1 , Nf+1 )
    Pnoise = get_tdiNSD( tditype , tdigen , tdiI , tdiJ , f )
    norm = np.sqrt( N / ( 2*stime ) * Pnoise )
    rez = np.random.standard_normal( ( Nf , ) )
    imz = np.random.standard_normal( ( Nf , ) )
    z = ( rez + 1j * imz ) / np.sqrt( 2 )
    ntilde_positivef = norm * z
    if N % 2 == 0 :
        ntilde = np.array( [0] + list( ntilde_positivef ) + [0] + list( np.flipud(np.conj(ntilde_positivef)) ) )
    else :
        ntilde = np.array( [0] + list( ntilde_positivef ) + list( np.flipud(np.conj(ntilde_positivef)) ) )
    n = np.real( sp.ifft( ntilde ) )
    return t , n


def get_tdiNoise_freq_domain_CovarMatrix( tditype , tdigen , whichtdis ,
                                          stime , duration , inittime , seed ) :
    """
    makes noise time-series given their covariance matrix
    INPUT:
    tditype --- type of TDI observable
    tdigen --- generation of TDI
    whichtdis --- this can either be '123' or 'AET'
    stime --- sampling time of the noise time-series
    duration --- duration of the noise time-series
    inittime --- initial time of the noise time-series
    seed --- seed for the random number generator
    OUTPUT:
    t --- time of the noise time-series
    n --- noise time-series
    """
    np.random.seed( seed )

    N = np.round( duration / stime )
    if N % 2 == 0 :
        Nf = int( N/2 - 1 )
    else :
        Nf = int( ( N-1 ) / 2 )
    df = 1 / (N*stime)
    f = df * np.arange( 1 , Nf + 1 )
    if whichtdis == '123' :
        pairs = [ ('1','1') , ('1','2') , ('1','3') ,
                  ('2','1') , ('2','2') , ('2','3') ,
                  ('3','1') , ('3','2') , ('3','3') ]
        cIJs = np.array( [ get_tdiNSD( tditype , tdigen , pair[0] , pair[1] , f ) for pair in pairs ] )
    elif whichtdis == 'AET' :
        pairs = [ ('A','A') , ('A','E') , ('A','T') ,
                  ('E','A') , ('E','E') , ('E','T') ,
                  ('T','A') , ('T','E') , ('T','T') ]
        cIJs = np.array( [ get_tdiNSD( tditype , tdigen , pair[0] , pair[1] , f ) for pair in pairs ] )
    comatrix = np.array( [ [ cIJs[0] , cIJs[1] , cIJs[2] ] , [ cIJs[3] , cIJs[4] , cIJs[5] ] , [ cIJs[6] , cIJs[7] , cIJs[8] ] ] )

    zs = np.array( [ ( np.random.standard_normal((Nf,)) + 1j * np.random.standard_normal((Nf,)) ) / np.sqrt(2)
                     for i in range( 3 ) ] )
    ntilde_p = np.zeros( ( 3 , Nf ) , dtype=complex )
    for k in range( Nf ) :
        C = comatrix[ :,:,k ]
        if not np.allclose( C , np.conj( np.transpose( C ) ) ) :
            print "Covariance matrix NOT Hermitian! Unphysical."        
        w , V = sp_linalg.eig( C )
        ntilde_p[ :,k ] =  np.sqrt( N / (2*stime) ) * np.dot( V , np.dot( np.sqrt( np.diag( w ) ) , zs[ :,k ] ) )
    
    zerofill = np.array( [ [0] , [0] , [0] ] )
    if N % 2 == 0 :
        ntilde = np.concatenate( ( zerofill , ntilde_p , zerofill , np.conj(np.fliplr(ntilde_p)) ) , axis = 1 )
    else :
        ntilde = np.concatenate( ( zerofill , ntilde_p , np.conj(np.fliplr(ntilde_p)) ) , axis = 1 )
    n = np.real( sp.ifft( ntilde , axis = 1 ) )
    t = inittime + stime * np.arange( N )
    return t , n
